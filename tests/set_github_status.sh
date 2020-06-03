# arguments:
# 1: GIT_COMMIT
# 2: GITHUB_TOKEN
# 3: RUN_DISPLAY_URL
# 4: status
echo "\"https://api.github.com/repos/Maarten-vd-Sande/snakemake-workflows/statuses/$1?access_token=$2\""
curl "https://api.github.com/repos/Maarten-vd-Sande/snakemake-workflows/statuses/${1}?access_token=${2}" \
-H "Content-Type: application/json" \
-X POST \
-d "{\"state\": \"$4\", \"description\": \"Jenkins\", \"target_url\": \"$3\"}"
