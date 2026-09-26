"""Create a numerical result without replacing an existing artifact."""
import json
from pathlib import Path


def write_result(path, result):
    """Reject invalid JSON and existing files before any result is replaced."""
    text = json.dumps(result, indent=2, allow_nan=False) + '\n'
    with Path(path).open('x') as stream:
        stream.write(text)
