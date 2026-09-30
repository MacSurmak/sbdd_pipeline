"""
File: cli/core/db.py
Description: SBDD Project State Database (LIMS). Single Source of Truth.
Uses SQLAlchemy for SQLite with future-proof compatibility for PostgreSQL.
"""

import json
from pathlib import Path
from sqlalchemy import (
    create_engine,
    Column,
    Integer,
    String,
    Boolean,
    ForeignKey,
    Text,
    Float,
)
from sqlalchemy.orm import declarative_base, sessionmaker, relationship

Base = declarative_base()


class Receptor(Base):
    """
    Stores macromolecular structures: raw crystals, AF3 models, MD frames, or fully prepared structures.
    """

    __tablename__ = "receptors"

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True, index=True, nullable=False)
    path = Column(String, nullable=False)
    status = Column(
        String, default="unprepared", nullable=False
    )  # unprepared, ready, error
    prep_mode = Column(String, default="full", nullable=False)  # full, light, ready
    parent_name = Column(String, ForeignKey("receptors.name"), nullable=True)
    meta_data = Column(Text, default="{}")  # JSON storage for cofactors, etc.

    parent = relationship("Receptor", remote_side=[name], backref="children")

    @property
    def meta_dict(self):
        return json.loads(self.meta_data)

    @meta_dict.setter
    def meta_dict(self, val):
        self.meta_data = json.dumps(val)


class Ligand(Base):
    """
    Stores ligand libraries or native reference ligands.
    """

    __tablename__ = "ligands"

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True, index=True, nullable=False)
    path = Column(String, nullable=False)
    status = Column(String, default="unprepared", nullable=False)
    is_native = Column(Boolean, default=False, nullable=False)


class Trajectory(Base):
    """
    Stores results of Molecular Dynamics simulations.
    Tracks both raw and cleaned (PBC whole, centered, fitted, solvent-stripped) trajectories.
    """

    __tablename__ = "trajectories"

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True, index=True, nullable=False)
    receptor_name = Column(String, ForeignKey("receptors.name"), nullable=False)
    path_xtc = Column(String, nullable=False)
    path_tpr = Column(String, nullable=False)
    path_clean_xtc = Column(String, nullable=True)
    path_ref_pdb = Column(String, nullable=True)
    length_ns = Column(Float, nullable=True)
    status = Column(String, default="raw", nullable=False)  # raw, cleaned, error

    receptor = relationship("Receptor")


class PoseSet(Base):
    """
    Stores massive screening/docking results (Library of poses).
    """

    __tablename__ = "pose_sets"

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True, index=True, nullable=False)
    path_sdf = Column(String, nullable=False)
    path_csv = Column(String, nullable=False)
    status = Column(String, default="ready", nullable=False)


class Complex(Base):
    """
    A specific Protein-Ligand hit chosen for downstream MD.
    """

    __tablename__ = "complexes"

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True, index=True, nullable=False)
    receptor_name = Column(String, ForeignKey("receptors.name"), nullable=False)
    ligand_title = Column(String, nullable=False)
    path_pose = Column(String, nullable=False)
    status = Column(String, default="ready", nullable=False)

    receptor = relationship("Receptor")


def init_db(db_path: Path = None):
    """
    Initializes the SQLite database. Creates tables if they don't exist.
    """
    if db_path is None:
        db_path = Path.cwd() / "sbdd_state.db"

    engine = create_engine(f"sqlite:///{db_path}")
    Base.metadata.create_all(engine)
    SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
    return SessionLocal()
