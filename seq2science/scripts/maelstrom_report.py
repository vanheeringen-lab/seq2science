import re
import base64
from contextlib import redirect_stdout, redirect_stderr

with open(str(snakemake.log), "w") as f:
   with redirect_stdout(f), redirect_stderr(f):
        png2base = lambda x: base64.b64encode(open(x, "rb").read()).decode('utf-8')

        with open(f"{snakemake.input[0]}/gimme.maelstrom.report.html", "r") as f:
            html = f.read()

        # search all png images and inject them into the report
        while match := re.search(
                "\<img\ src\=(\"logos\/(.+)\.png\")\ style\=\"height\:100\%\;width\:100\%\;object-fit\:contain\;\"\/>", html):
            png_loc = f"{snakemake.input[0]}/logos/{match.group(2)}.png"
            png_as_base64 = png2base(png_loc)
            inject_png = f"data:image/png;base64,{png_as_base64}\""

            start, stop = match.span()
            start += len("<img src =")
            stop = start + len(match.group(1))

            html = html[:start] + inject_png + html[stop:]

        # make the logo rows slightly larger
        html = re.sub("<div style=\"height:30px;object-fit:contain;\">",
                      "<div style=\"height:50px;object-fit:contain;\">",
                      html)
        html = re.sub("<th class=\"col_heading level0 col1\" >logo</th>",
                      "<th class=\"col_heading level0 col1\" >motif information</th>",
                      html)


        with open(snakemake.output, "w") as f:
            f.write(html)
