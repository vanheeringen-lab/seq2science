import re
from setuptools import setup
from distutils.util import convert_path
import ast
import toml

project = toml.load("pyproject.toml")["project"]

# read the readme as long description
with open("README.md") as f:
    project["long_description"] = f.read()

with open(convert_path('seq2science/__init__.py')) as ver_file:
    match = next(re.finditer('__version__ = "(.*)"', ver_file.read(), re.MULTILINE))
    project["version"] = match.group(1)
project["long_description_content_type"] = "text/markdown"
project["package_data"] = ast.literal_eval(project["package_data"])

setup(**project)
