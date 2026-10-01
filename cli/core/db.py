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
    Stores individual chemical molecules (2D/3D structure without receptor positioning).
    """

    __tablename__ = "ligands"

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True, index=True, nullable=False)
    path = Column(String, nullable=False)
    smiles = Column(Text, nullable=True)
    status = Column(String, default="unprepared", nullable=False)  # unprepared, ready
    meta_data = Column(Text, default="{}")


class Pose(Base):
    """
    A specific 3D spatial conformation of a Ligand positioned in a Receptor coordinate frame.
    """

    __tablename__ = "poses"

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True, index=True, nullable=False)
    ligand_name = Column(String, ForeignKey("ligands.name"), nullable=False)
    receptor_name = Column(String, ForeignKey("receptors.name"), nullable=False)
    origin = Column(String, default="docking", nullable=False)  # crystal, docking, manual
    path = Column(String, nullable=False)
    status = Column(String, default="unprepared", nullable=False)  # unprepared, ready
    scores = Column(Text, default="{}")

    ligand = relationship("Ligand")
    receptor = relationship("Receptor")


class System(Base):
    """
    The physical thermodynamic unit ready for MD or Scoring: Receptor + optional Pose (None = Apo).
    Supports hierarchical conformational ensembles via parent_system_name.
    """

    __tablename__ = "systems"

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True, index=True, nullable=False)
    receptor_name = Column(String, ForeignKey("receptors.name"), nullable=False)
    pose_name = Column(String, ForeignKey("poses.name"), nullable=True)
    parent_system_name = Column(String, ForeignKey("systems.name"), nullable=True)
    origin = Column(String, default="crystal", nullable=False)  # crystal, docking, md_cluster
    status = Column(String, default="unprepared", nullable=False)  # unprepared, ready
    meta_data = Column(Text, default="{}")

    receptor = relationship("Receptor")
    pose = relationship("Pose")
    parent_system = relationship(
        "System", remote_side=[name], backref="conformations"
    )


class Trajectory(Base):
    """
    Stores Molecular Dynamics simulations linked to a physical System.
    """

    __tablename__ = "trajectories"

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True, index=True, nullable=False)
    system_name = Column(String, ForeignKey("systems.name"), nullable=False)
    replicate = Column(Integer, default=1, nullable=False)
    path_xtc = Column(String, nullable=False)
    path_tpr = Column(String, nullable=False)
    path_clean_xtc = Column(String, nullable=True)
    path_ref_pdb = Column(String, nullable=True)
    length_ns = Column(Float, nullable=True)
    status = Column(String, default="raw", nullable=False)  # raw, cleaned, error

    system = relationship("System")


class PoseSet(Base):
    """
    Grouping table for virtual screening results.
    """

    __tablename__ = "pose_sets"

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True, index=True, nullable=False)
    path_sdf = Column(String, nullable=False)
    path_csv = Column(String, nullable=False)
    status = Column(String, default="ready", nullable=False)


def init_db(db_path: Path = None):
    """
    Initializes the SQLite database. Creates tables if they don't exist.
    """
    if db_path is None:
        db_path = Path.cwd() / "sbdd_state.db"

    engine = create_engine(f"sqlite:///{db_path}")
    Base.metadata.create_all(engine)
    SessionLocal = sessionmaker(
        autocommit=False, autoflush=False, expire_on_commit=False, bind=engine
    )
    return SessionLocal()
