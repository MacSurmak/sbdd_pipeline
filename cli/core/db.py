"""
File: cli/core/db.py
Description: Clean 3NF SBDD Project Relational Database Schema.
Enforces surrogate integer primary keys, foreign key constraints, and state lineage.
Zero string-based foreign keys. Zero suffix renaming hacks.
"""

import json
from pathlib import Path
from datetime import datetime
from sqlalchemy import (
    create_engine,
    event,
    Column,
    Integer,
    String,
    ForeignKey,
    Text,
    Float,
    DateTime,
)
from sqlalchemy.engine import Engine
from sqlalchemy.orm import declarative_base, sessionmaker, relationship

Base = declarative_base()


@event.listens_for(Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    """Enforce foreign key constraints in SQLite."""
    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA foreign_keys=ON")
    cursor.close()


class Target(Base):
    """
    Biological macromolecular target / project folder (e.g. 'RyR1', 'TSPO', 'EGFR').
    Contains no coordinate files; acts strictly as an abstract metadata container.
    """

    __tablename__ = "targets"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String, unique=True, index=True, nullable=False)
    description = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    frames = relationship(
        "ReceptorFrame", back_populates="target", cascade="all, delete-orphan"
    )
    systems = relationship(
        "System", back_populates="target", cascade="all, delete-orphan"
    )


class ReceptorFrame(Base):
    """
    Specific 3D structure of a Target (crystal model, NMR model, MD centroid, AlphaFold).
    All frames are peer entities under their parent Target.
    """

    __tablename__ = "receptor_frames"

    id = Column(Integer, primary_key=True, autoincrement=True)
    target_id = Column(
        Integer, ForeignKey("targets.id", ondelete="CASCADE"), nullable=False
    )
    name = Column(String, unique=True, index=True, nullable=False)
    source_type = Column(
        String, default="crystal", nullable=False
    )  # crystal, nmr, md_centroid, af3, docking_receptor
    state = Column(String, default="raw", nullable=False)  # raw, ready
    parent_frame_id = Column(
        Integer, ForeignKey("receptor_frames.id", ondelete="SET NULL"), nullable=True
    )
    path = Column(String, nullable=False)
    cofactors = Column(Text, default="[]", nullable=False)  # JSON list of cofactors
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    target = relationship("Target", back_populates="frames")
    parent_frame = relationship(
        "ReceptorFrame", remote_side=[id], backref="derived_frames"
    )
    poses = relationship(
        "Pose", back_populates="receptor_frame", cascade="all, delete-orphan"
    )
    systems = relationship(
        "System", back_populates="receptor_frame", cascade="all, delete-orphan"
    )

    @property
    def cofactors_list(self) -> list:
        try:
            return json.loads(self.cofactors) if self.cofactors else []
        except Exception:
            return []

    @cofactors_list.setter
    def cofactors_list(self, val: list):
        self.cofactors = json.dumps(val if val is not None else [])


class Ligand(Base):
    """
    Chemical substance (2D topology).
    Deduplicated across the entire database by Standard InChIKey.
    """

    __tablename__ = "ligands"

    id = Column(Integer, primary_key=True, autoincrement=True)
    inchikey = Column(String(27), unique=True, index=True, nullable=False)
    smiles = Column(Text, nullable=False)
    name = Column(String, nullable=False)
    batch_id = Column(String, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    poses = relationship("Pose", back_populates="ligand", cascade="all, delete-orphan")


class Pose(Base):
    """
    Atomic 3D conformation of a Ligand positioned in a ReceptorFrame pocket.
    1 Pose = 1 File on disk (.sdf or .mae).
    """

    __tablename__ = "poses"

    id = Column(Integer, primary_key=True, autoincrement=True)
    ligand_id = Column(
        Integer, ForeignKey("ligands.id", ondelete="CASCADE"), nullable=False
    )
    receptor_frame_id = Column(
        Integer, ForeignKey("receptor_frames.id", ondelete="CASCADE"), nullable=False
    )
    name = Column(String, unique=True, index=True, nullable=False)
    origin = Column(
        String, default="crystal", nullable=False
    )  # crystal, docking, md_cluster, manual
    state = Column(String, default="raw", nullable=False)  # raw, ready
    parent_pose_id = Column(
        Integer, ForeignKey("poses.id", ondelete="SET NULL"), nullable=True
    )
    path = Column(String, nullable=False)
    scores = Column(Text, default="{}", nullable=False)  # JSON dict of scoring metrics
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    ligand = relationship("Ligand", back_populates="poses")
    receptor_frame = relationship("ReceptorFrame", back_populates="poses")
    parent_pose = relationship("Pose", remote_side=[id], backref="derived_poses")
    systems = relationship("System", back_populates="pose")

    @property
    def scores_dict(self) -> dict:
        try:
            return json.loads(self.scores) if self.scores else {}
        except Exception:
            return {}

    @scores_dict.setter
    def scores_dict(self, val: dict):
        self.scores = json.dumps(val if val is not None else {})


class System(Base):
    """
    Thermodynamic physical unit for simulation or evaluation:
    ReceptorFrame + optional Pose (NULL = Apo, Not NULL = Holo).
    """

    __tablename__ = "systems"

    id = Column(Integer, primary_key=True, autoincrement=True)
    target_id = Column(
        Integer, ForeignKey("targets.id", ondelete="CASCADE"), nullable=False
    )
    receptor_frame_id = Column(
        Integer, ForeignKey("receptor_frames.id", ondelete="CASCADE"), nullable=False
    )
    pose_id = Column(
        Integer, ForeignKey("poses.id", ondelete="SET NULL"), nullable=True
    )
    name = Column(String, unique=True, index=True, nullable=False)
    state = Column(String, default="raw", nullable=False)  # raw, ready
    origin = Column(
        String, default="crystal", nullable=False
    )  # crystal, docking, md_cluster
    environment = Column(
        Text, default="{}", nullable=False
    )  # JSON membrane/solvent specs
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    target = relationship("Target", back_populates="systems")
    receptor_frame = relationship("ReceptorFrame", back_populates="systems")
    pose = relationship("Pose", back_populates="systems")
    trajectories = relationship(
        "Trajectory", back_populates="system", cascade="all, delete-orphan"
    )

    @property
    def environment_dict(self) -> dict:
        try:
            return json.loads(self.environment) if self.environment else {}
        except Exception:
            return {}

    @environment_dict.setter
    def environment_dict(self, val: dict):
        self.environment = json.dumps(val if val is not None else {})


class Trajectory(Base):
    """
    Molecular Dynamics simulation replicate linked to a physical System.
    Raw chunks live in raw_dir; cleaned trajectory without solvent lives in path_clean_xtc.
    """

    __tablename__ = "trajectories"

    id = Column(Integer, primary_key=True, autoincrement=True)
    system_id = Column(
        Integer, ForeignKey("systems.id", ondelete="CASCADE"), nullable=False
    )
    name = Column(String, unique=True, index=True, nullable=False)
    replicate = Column(Integer, default=1, nullable=False)
    length_ns = Column(Float, default=0.0, nullable=False)
    raw_dir = Column(String, nullable=False)
    path_tpr = Column(String, nullable=False)
    path_clean_xtc = Column(String, nullable=True)
    path_ref_pdb = Column(String, nullable=True)
    status = Column(String, default="raw", nullable=False)  # raw, cleaned, error
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    system = relationship("System", back_populates="trajectories")


@event.listens_for(ReceptorFrame.path, "set")
def forbid_scratch_in_receptor_frame_path(target, value, oldvalue, initiator):
    if value and ("/work/" in str(value) or ".nextflow" in str(value)):
        raise ValueError(
            f"FORBIDDEN: Attempted to write ephemeral scratch path into ReceptorFrame.path: {value}. "
            f"All database records must strictly reference permanent project deliverables."
        )


@event.listens_for(Pose.path, "set")
def forbid_scratch_in_pose_path(target, value, oldvalue, initiator):
    if value and ("/work/" in str(value) or ".nextflow" in str(value)):
        raise ValueError(
            f"FORBIDDEN: Attempted to write ephemeral scratch path into Pose.path: {value}. "
            f"All database records must strictly reference permanent project deliverables."
        )


def init_db(db_path: Path = None):
    """
    Initializes SQLite connection and creates all 3NF tables if not existing.
    """
    if db_path is None:
        db_path = Path.cwd() / "sbdd_state.db"

    engine = create_engine(f"sqlite:///{db_path}")
    Base.metadata.create_all(engine)
    SessionLocal = sessionmaker(
        autocommit=False, autoflush=False, expire_on_commit=False, bind=engine
    )
    return SessionLocal()
