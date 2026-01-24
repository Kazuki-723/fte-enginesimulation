import json
import re
from pathlib import Path

def load_jsonc(jsonc_path):
    """JSONC ファイルを読み込み、コメントを除去して dict を返す"""
    jsonc_path = Path(jsonc_path)
    if not jsonc_path.exists():
        raise FileNotFoundError(f"cannot found jsonc file: {jsonc_path}")

    text = jsonc_path.read_text(encoding="utf-8")

    # // コメント削除
    text = re.sub(r"//.*", "", text)
    # /* ... */ コメント削除
    text = re.sub(r"/\*[\s\S]*?\*/", "", text)

    return json.loads(text)


if __name__ == "__main__":
    print("input jsonc file name, 'example.jsonc' :")
    filename = input("> ").strip()

    try:
        config = load_jsonc(filename)
    except Exception as e:
        print(f"loading error: {e}")
        exit(1)

    # --- JSONC の値を変数に割り当て ---
    a = config.get("a")
    b = config.get("b")

    print("input values:")
    print(f"a = {a}")
    print(f"b = {b}")

    # 計算例
    if a is not None and b is not None:
        print(f"a + b = {a + b}")