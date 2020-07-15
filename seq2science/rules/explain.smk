import re
import os.path
import yaml


def explain_rule(string):
    """
    Parse a message. Depending on whether explain_rule is true, we do nothing
    or we parse.
    """
    if config.get("explain_rule") != True:
        return None
    else:
        # clean our explanation
        string = string.replace("\n", "")
        string = " ".join(string.split())

        # find our environment
        env_dir = os.path.normpath(os.path.join(config['rule_dir'], "..", "envs"))

        parser = re.compile("@([^[]*)\[([^]]*)\]")
        while len(parser.findall(string)):
            match = next(parser.finditer(string))

            # parse our environment
            yaml_file = f"{env_dir}/{match.group(1)}.yaml"
            tool = match.group(2)
            with open(yaml_file, 'r') as stream:
                env = yaml.safe_load(stream)

            for dependency in env["dependencies"]:
                if tool in dependency:
                    version = dependency[dependency.find("=") + 1:]
                    break
            else:
                continue

            # replace the placeholder with the actual version
            string = string[:match.span()[0]] + version + string[match.span()[1]:]
        return string
