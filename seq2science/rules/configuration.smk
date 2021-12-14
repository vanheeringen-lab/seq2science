"""
The seq2science configuration/preprocessing is split into four parts:
* generic: all logic not related to any specific workflows
* workflows: all logic related to specific workflows
* explain: all logic necessary to make an explanation of what has/will be done
* logging: all logic related to logging to stdout/file
"""


include: "../rules/configuration_generic.smk"
include: "../rules/configuration_workflows.smk"
include: "../rules/configuration_explain.smk"
include: "../rules/configuration_logging.smk"
