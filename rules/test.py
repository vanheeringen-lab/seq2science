import conda.cli.python_api

stdout, stderr, return_code = conda.cli.python_api.run_command('search', 'r-base', channel='r')

versions = []
for line in stdout.split('\n'):
    if not line == '':
        output = [sub for sub in line.split(' ') if sub != '']

        if output[-1] == 'conda-forge':
            versions.append(output)

for v in versions:
    print(v)
# print(versions)