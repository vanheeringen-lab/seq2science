from setuptools import setup, find_packages
import ast
import toml

project = toml.load("pyproject.toml")["project"]

# read the readme as long description
with open("README.md") as f:
    project["long_description"] = f.read()

project["long_description_content_type"] = "text/markdown"
project["package_data"] = ast.literal_eval(project["package_data"])
import sys
sys.stdout.write("\n\n\n\n")
sys.stdout.write(str(project))
sys.stdout.write("\n\n\n\n")
# find packages to install
# project["packages"] = find_packages(where='./*', exclude=())
# # print()
# assert False, project["packages"]

setup(**project)
