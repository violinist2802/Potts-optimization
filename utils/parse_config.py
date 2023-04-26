import yaml
def parse_config(path: str, scenario: str) -> dict:
    """
    Parse yaml file

    Args:
        path: str, path to config

    Returns:
        cfg: dict, dictionary with new params
    """
    with open(path) as f:
        data= yaml.load(f, Loader=yaml.FullLoader)
        cfg=data['neonatal_rat'][scenario]

    return cfg
