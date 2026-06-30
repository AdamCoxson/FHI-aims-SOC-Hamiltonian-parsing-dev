# Physical Consistency Unit Tests

This directory contains `pytest` unit tests designed to verify the physical consistency and symmetries of the noncollinear Spin-Orbit Coupling (SOC) Hamiltonians and Overlap matrices parsed from FHI-aims.

## Prerequisites

You need `pytest` installed in your active environment to run these scripts:
```bash
conda activate fhi-deephdock
```
*(Note: `pytest` should already be installed in this environment from earlier steps).*

## Running the Tests

To run the full test suite, simply point `pytest` to this directory:
```bash
pytest ./
```

### Running Specific Tests
You can also execute a specific test file directly:
```bash
pytest test_overlap.py
```

## Configuring the Data Directory

By default, the tests are configured (via `conftest.py`) to look for the intermediate `.h5` files in the following relative path:
`../../2D_materials_data/MoS2_unitcell/LAO_processed_data/deeph_way/soc_intermediate`

### Overriding the Data Directory on the Fly

If you want to test matrices from a different parser run (for example, your `alternative_way` output), you don't need to edit any Python code. Instead, use the custom `--data-dir` command-line argument:

```bash
pytest ./ --data-dir=../../2D_materials_data/MoS2_unitcell/LAO_processed_data/alternative_way/soc_intermediate
```

### Changing the Default Directory

If you ever want to permanently change the default fallback directory, just open `conftest.py` and modify the `default_path` variable on Line 5:

```python
def pytest_addoption(parser):
    # Edit the path below to change the default fallback directory
    default_path = Path(__file__).parent.parent.parent / "2D_materials_data" / "MoS2_unitcell" / "LAO_processed_data" / "alternative_way" / "soc_intermediate"
    ...
```
