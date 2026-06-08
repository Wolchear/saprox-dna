from typing import TypedDict

class Config(TypedDict):
    base_root: str
    dirs: dict[str, str]

def get_path(config: Config, key: str) -> str:
    return f"{config['base_root']}/{config['dirs'][key]}"