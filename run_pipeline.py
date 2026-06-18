from __future__ import annotations

import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))

main = importlib.import_module('modelyml.driver').main


if __name__ == '__main__':
    raise SystemExit(main())
