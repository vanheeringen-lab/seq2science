import re
import os
import ast
import toml
import shutil
from setuptools import setup
from distutils.util import convert_path
from setuptools.command.install import install


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        # move the requirements.yaml to the seq2science environments
        env_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "requirements.yaml")
        build_dir = os.path.join(self.build_lib, "seq2science", "envs", "seq2science_requirements.yaml")
        os.makedirs(os.path.dirname(build_dir), exist_ok=True)
        shutil.copy(env_file, build_dir)

        # install as normally
        install.run(self)


project = toml.load("pyproject.toml")["project"]

# read the readme as long description
with open("README.md") as f:
    project["long_description"] = f.read()

with open(convert_path('seq2science/__init__.py')) as ver_file:
    match = next(re.finditer('__version__ = "(.*)"', ver_file.read(), re.MULTILINE))
    project["version"] = match.group(1)
project["long_description_content_type"] = "text/markdown"
project["package_data"] = ast.literal_eval(project["package_data"])
project["cmdclass"] = {
    "install": PostInstallCommand,
}

setup(**project)
