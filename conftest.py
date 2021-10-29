import pytest
from pathlib import Path


@pytest.fixture
def refs(request):
    return Path(request.module.__file__).parent / 'refs'
