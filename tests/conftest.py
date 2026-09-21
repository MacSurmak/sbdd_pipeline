import pytest
from pathlib import Path

SAMPLE_SDF = """MOL1
     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
> <CNNscore>
0.8500

> <CNNaffinity>
7.200

> <minimizedAffinity>
-9.500

$$$$
MOL2
     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
> <CNNscore>
0.6200

> <CNNaffinity>
6.100

> <minimizedAffinity>
-8.100

$$$$
MOL3
     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
> <CNNscore>
0.9100

> <CNNaffinity>
8.400

> <minimizedAffinity>
-11.200

$$$$
"""

SAMPLE_DATAPAR1 = "TSPO_orchestra;20.7; 0.9;  28;   2; -20.2;   0.0;\n"

@pytest.fixture
def sample_sdf(tmp_path):
    sdf_file = tmp_path / "test_lib.sdf"
    sdf_file.write_text(SAMPLE_SDF)
    return sdf_file

@pytest.fixture
def sample_datapar1(tmp_path):
    dat_file = tmp_path / "datapar1"
    dat_file.write_text(SAMPLE_DATAPAR1)
    return dat_file

