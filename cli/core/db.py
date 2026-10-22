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
    Table,
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


system_cofactors = Table(
    "system_cofactors",
    Base.metadata,
    Column(
        "system_id",
        Integer,
        ForeignKey("systems.id", ondelete="CASCADE"),
        primary_key=True,
    ),
    Column(
        "pose_id", Integer, ForeignKey("poses.id", ondelete="CASCADE"), primary_key=True
    ),
)


class Target(Base):
    """
    Biological macromolecular target / project folder (e.g. 'RyR1', 'TSPO', 'EGFR').
    Contains no coordinate files; acts strictly as an abstract metadata container.
    """

    __tablename__ = "targets"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String, unique=True, index=True, nullable=False)
    reference_sequence = Column(Text, nullable=True)  # SSOT Target FASTA sequence
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
    Stores explicit, validated paths to both Maestro MAE and PDB deliverables.
    """

    __tablename__ = "receptor_frames"

    id = Column(Integer, primary_key=True, autoincrement=True)
    target_id = Column(
        Integer, ForeignKey("targets.id", ondelete="CASCADE"), nullable=False
    )
    name = Column(String, unique=True, index=True, nullable=False)
    source_type = Column(
        String, default="crystal", nullable=False
    )  # crystal, nmr, md_centroid, docking_receptor, crystal_healed
    state = Column(String, default="raw", nullable=False)  # broken, raw, ready
    parent_frame_id = Column(
        Integer, ForeignKey("receptor_frames.id", ondelete="SET NULL"), nullable=True
    )
    source_trajectory_id = Column(
        Integer, ForeignKey("trajectories.id", ondelete="SET NULL"), nullable=True
    )
    feature_group = Column(String(64), nullable=True)  # Selection group for coordinates (e.g. Pocket, Backbone, Ligand)
    align_group = Column(String(64), nullable=True)    # Basis group for superimposition (e.g. Backbone, C-alpha)
    dim_reduction_method = Column(String(16), nullable=True)  # 'pca' or 'tica'
    cluster_rank = Column(Integer, nullable=True)      # Medoid population ranking (1, 2, 3...)
    population = Column(Float, nullable=True)          # Exact Boltzmann ensemble fraction (0.0 to 1.0)
    compactness_spread = Column(Float, nullable=True)  # Intra-cluster coordinate spread (Å)
    path = Column(String, nullable=False)  # Canonical primary path
    path_mae = Column(
        String, nullable=True
    )  # Validated MAE deliverable for Glide/Prime
    path_pdb = Column(
        String, nullable=True
    )  # Validated PDB deliverable for GNINA/tLEaP
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    target = relationship("Target", back_populates="frames")
    source_trajectory = relationship("Trajectory", foreign_keys=[source_trajectory_id])
    parent_frame = relationship(
        "ReceptorFrame", remote_side=[id], backref="derived_frames"
    )
    poses = relationship(
        "Pose", back_populates="receptor_frame", cascade="all, delete-orphan"
    )
    systems = relationship(
        "System", back_populates="receptor_frame", cascade="all, delete-orphan"
    )


class Substance(Base):
    """
    Abstract Chemical Substance (Parent 2D compound).
    Deduplicated across the entire database by Neutral Standard InChIKey.
    Acts strictly as a catalogue-level container for all chemical microstates.
    """

    __tablename__ = "substances"

    id = Column(Integer, primary_key=True, autoincrement=True)
    inchikey = Column(String(27), unique=True, index=True, nullable=False)
    smiles = Column(Text, nullable=False)
    name = Column(String, nullable=False)
    res_name = Column(String(8), nullable=False)
    substance_type = Column(
        String, default="ligand", nullable=False
    )  # 'ligand', 'cofactor', 'ion'
    state = Column(String, default="raw", nullable=False)  # 'raw', 'ready'
    batch_id = Column(String, nullable=True)
    path_mae = Column(String, nullable=True)  # Prepared .maegz for LigPrep
    path_sdf = Column(String, nullable=True)  # Prepared .sdf
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    states = relationship(
        "SubstanceState", back_populates="substance", cascade="all, delete-orphan"
    )


class SubstanceState(Base):
    """
    Specific Chemical Microstate (Protomer, Tautomer, Stereoisomer, Ionization state).
    Keyed deterministically by state_hash (SHA-256 of RDKit Isomeric Canonical SMILES).
    Acts as the single quantum invariant for GAFF2 force field parameters and AM1-BCC charges.
    """

    __tablename__ = "substance_states"

    id = Column(Integer, primary_key=True, autoincrement=True)
    substance_id = Column(
        Integer, ForeignKey("substances.id", ondelete="CASCADE"), nullable=False
    )
    state_hash = Column(String(32), index=True, nullable=False)
    isomeric_smiles = Column(Text, nullable=False)
    formal_charge = Column(
        Integer, nullable=False
    )  # Exact net charge for -nc in Antechamber
    state = Column(String, default="raw", nullable=False)  # 'raw', 'ready'
    path_sdf = Column(
        String, nullable=True
    )  # 3D reference conformer for this specific state
    path_mol2 = Column(String, nullable=True)  # GAFF2 AM1-BCC topology
    path_frcmod = Column(String, nullable=True)  # GAFF2 force field parameters
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    substance = relationship("Substance", back_populates="states")
    poses = relationship(
        "Pose", back_populates="substance_state", cascade="all, delete-orphan"
    )


class Pose(Base):
    """
    Spatial 3D conformation of a specific SubstanceState in a ReceptorFrame pocket.
    Contains Cartesian coordinates (x, y, z) while delegating chemical topology to SubstanceState.
    """

    __tablename__ = "poses"

    id = Column(Integer, primary_key=True, autoincrement=True)
    substance_state_id = Column(
        Integer, ForeignKey("substance_states.id", ondelete="CASCADE"), nullable=False
    )
    receptor_frame_id = Column(
        Integer, ForeignKey("receptor_frames.id", ondelete="CASCADE"), nullable=False
    )
    name = Column(String, unique=True, index=True, nullable=False)
    res_name = Column(String(8), nullable=False)
    origin = Column(
        String, default="crystal", nullable=False
    )  # crystal, docking_glide, docking_gnina, md_cluster
    state = Column(String, default="ready", nullable=False)  # raw, ready
    parent_pose_id = Column(
        Integer, ForeignKey("poses.id", ondelete="SET NULL"), nullable=True
    )
    path = Column(String, nullable=False)  # Canonical primary path
    path_sdf = Column(String, nullable=True)  # Validated SDF deliverable
    path_mae = Column(String, nullable=True)  # Validated MAE deliverable
    scores = Column(Text, default="{}", nullable=False)  # JSON dict of scoring metrics
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    substance_state = relationship("SubstanceState", back_populates="poses")
    receptor_frame = relationship("ReceptorFrame", back_populates="poses")
    parent_pose = relationship("Pose", remote_side=[id], backref="derived_poses")
    systems = relationship("System", back_populates="pose")

    @property
    def substance(self):
        return self.substance_state.substance if self.substance_state else None

    @property
    def formal_charge(self) -> int:
        return self.substance_state.formal_charge if self.substance_state else 0

    @property
    def path_mol2(self) -> str | None:
        return self.substance_state.path_mol2 if self.substance_state else None

    @property
    def path_frcmod(self) -> str | None:
        return self.substance_state.path_frcmod if self.substance_state else None

    @property
    def scores_dict(self) -> dict:
        if not self.scores:
            return {}
        return json.loads(self.scores)

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
    path_gro = Column(String, nullable=True)  # Prepared canonical coordinate file
    path_top = Column(String, nullable=True)  # Prepared canonical topology file
    path_ndx = Column(String, nullable=True)  # Prepared canonical index file
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    target = relationship("Target", back_populates="systems")
    receptor_frame = relationship("ReceptorFrame", back_populates="systems")
    pose = relationship("Pose", back_populates="systems")
    cofactors = relationship(
        "Pose", secondary=system_cofactors, backref="cofactor_systems"
    )
    trajectories = relationship(
        "Trajectory", back_populates="system", cascade="all, delete-orphan"
    )

    @property
    def environment_dict(self) -> dict:
        if not self.environment:
            return {}
        return json.loads(self.environment)

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


def _validate_non_scratch_path(target, value, oldvalue, initiator):
    """
    Absolute barrier between ephemeral Nextflow/Slurm scratch and permanent DB records.
    Triggers fatal exception on any attempt to persist intermediate work directories.
    """
    if value:
        val_str = str(value)
        forbidden_tokens = ["/work/", ".nextflow", "/tmp/", "slurm_tmp", "nfs_scratch"]
        for token in forbidden_tokens:
            if token in val_str:
                col_name = initiator.key if hasattr(initiator, "key") else "path"
                cls_name = target.__class__.__name__
                raise ValueError(
                    f"FATAL DB INTEGRITY VIOLATION: Ephemeral scratch path detected in {cls_name}.{col_name}: '{val_str}'! "
                    f"Only validated, permanent project deliverables may be recorded in SQLite 3NF."
                )


# Attach strict validator to ALL file-referencing attributes across the entire schema
for _model_attr in [
    ReceptorFrame.path,
    ReceptorFrame.path_mae,
    ReceptorFrame.path_pdb,
    Substance.path_mae,
    Substance.path_sdf,
    SubstanceState.path_sdf,
    SubstanceState.path_mol2,
    SubstanceState.path_frcmod,
    Pose.path,
    Pose.path_sdf,
    Pose.path_mae,
    System.path_gro,
    System.path_top,
    System.path_ndx,
    Trajectory.raw_dir,
    Trajectory.path_tpr,
    Trajectory.path_clean_xtc,
    Trajectory.path_ref_pdb,
]:
    event.listen(_model_attr, "set", _validate_non_scratch_path)


FRAMEWORK_DIR = Path(__file__).resolve().parent.parent.parent


def init_db(db_path: Path = None, create_if_missing: bool = False):
    """
    Initializes SQLite connection to the 3NF database.
    Strictly forbids database operations inside the framework source tree (/shared/sbdd_pipeline).
    Enforces Fail-Fast if the database does not exist on query/execution commands.
    """
    if db_path is None:
        db_path = Path.cwd() / "sbdd_state.db"
    else:
        db_path = Path(db_path).resolve()

    # 1. Absolute barrier: Framework source repository is an immutable monolith
    if db_path.parent.resolve() == FRAMEWORK_DIR.resolve():
        raise PermissionError(
            f"FATAL: Attempted to create or access database inside pipeline codebase: {db_path}!\n"
            f"The framework repository at {FRAMEWORK_DIR} is an immutable monolith. "
            f"All SBDD projects must be initialized and executed in dedicated project directories (e.g. /shared/<project_name>)."
        )

    # 2. Fail-Fast: Refuse to silently spawn ghost databases on read/run commands
    if not create_if_missing and not db_path.exists():
        raise FileNotFoundError(
            f"FATAL: Project database not found at {db_path}!\n"
            f"Current working directory '{Path.cwd()}' is not an initialized SBDD project. "
            f"Run 'sbdd init' to initialize a project workspace first."
        )

    engine = create_engine(f"sqlite:///{db_path}")
    if create_if_missing:
        Base.metadata.create_all(engine)

    SessionLocal = sessionmaker(
        autocommit=False, autoflush=False, expire_on_commit=False, bind=engine
    )
    return SessionLocal()
