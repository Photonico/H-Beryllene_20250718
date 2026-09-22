import pytest


@pytest.fixture(
    params=["force_and_stress", "local_potential", "structure", "occupancies"]
)
def function_name(request):
    return request.param
