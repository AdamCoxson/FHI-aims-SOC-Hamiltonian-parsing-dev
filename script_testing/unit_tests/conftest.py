import pytest
from pathlib import Path

def pytest_addoption(parser):
    default_path = Path(__file__).parent.parent.parent / "2D_materials_data" / "MoS2_unitcell" / "LAO_processed_data" / "deeph_way" / "soc_intermediate"
    parser.addoption(
        "--data-dir", 
        action="store", 
        default=str(default_path), 
        help="Path to SOC intermediate h5 directory"
    )

@pytest.fixture
def data_dir(request):
    return Path(request.config.getoption("--data-dir"))
