"""Execute documentation examples and compare the printed output."""

import contextlib
import io
import re
import sys
import types
from pathlib import Path


def main():
    count = 0
    for path in sorted(Path("docs").rglob("*.md")):
        module = types.ModuleType("docs_example_" + path.stem.replace("-", "_"))
        sys.modules[module.__name__] = module
        blocks = list(re.finditer(r"```(python|text)\n(.*?)\n```", path.read_text(), re.S))
        for index, match in enumerate(blocks):
            if match[1] != "python" or index + 1 == len(blocks) or blocks[index + 1][1] != "text":
                continue
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                exec(compile(match[2], str(path), "exec"), module.__dict__)
            if index + 1 < len(blocks) and blocks[index + 1][1] == "text":
                expected = blocks[index + 1][2].rstrip()
                actual = output.getvalue().rstrip()
                if actual != expected:
                    raise AssertionError(f"{path}: {actual!r} != {expected!r}")
            count += 1
    print(f"Executed {count} documentation examples.")


if __name__ == "__main__":
    main()
