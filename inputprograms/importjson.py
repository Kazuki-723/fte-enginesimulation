import json
import re
from pathlib import Path

class JsoncLoader:
    """JSONC ファイルを読み込み、辞書として返すクラス"""

    def __init__(self, path: str):
        self.path = Path(path)

    def load(self) -> dict:
        """JSONC を読み込み、コメントを除去して dict を返す"""
        if not self.path.exists():
            raise FileNotFoundError(f"cannot found json file: {self.path}")

        text = self.path.read_text(encoding="utf-8")

        # // コメント削除
        text = re.sub(r"//.*", "", text)
        # /* ... */ コメント削除
        text = re.sub(r"/\*[\s\S]*?\*/", "", text)

        return json.loads(text)