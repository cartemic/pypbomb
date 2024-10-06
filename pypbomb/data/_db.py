import sqlite3
from pathlib import Path


def connection() -> sqlite3.Connection:
    return sqlite3.connect(Path(__file__).parent / "data.sqlite")
