import yaml
import os.path
from .compass import *  # noqa: F403
from .positions import *  # noqa: F403

datadir = os.path.expanduser(
    "~/work/COMPASS/Reduced/Imager/"
)  # assume that we are on ERDA
erda = True
try:
    with open("config.yaml") as f:
        config = yaml.safe_load(f)
    datadir = config["datadir"]
    erda = config["erda"]
except (FileNotFoundError, KeyError):
    pass
