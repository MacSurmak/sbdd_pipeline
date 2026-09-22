#!/usr/bin/env python3
import os
import sys
import re
import shutil
import argparse
import subprocess
from pathlib import Path

# --- НАСТРОЙКИ ИГНОРИРОВАНИЯ ДЛЯ ДАМПА ---
DEFAULT_IGNORE_DIRS = {
    '.git', '__pycache__', '.venv', 'venv', 'env', 
    'node_modules', '.idea', '.vscode', 'dist', 'build', 
    '.pytest_cache', '.mypy_cache', 'ppm3_server_code', '.nextflow',
    '.patch_backups'
}

DEFAULT_IGNORE_FILES = {
    '.DS_Store', 'Thumbs.db', 'devtool.py', 'dump.py', 'project_dump.txt'
}

TEXT_EXTENSIONS = {
    '.py', '.js', '.jsx', '.ts', '.tsx', '.html', '.css', '.scss', '.sass',
    '.json', '.yaml', '.yml', '.toml', '.xml', '.md', '.rst', '.txt',
    '.c', '.cpp', '.h', '.hpp', '.cs', '.java', '.go', '.rs', '.php',
    '.rb', '.sh', '.bat', '.ps1', '.sql', '.dockerfile', '.env.example',
    '.gitignore', '.ini', '.cfg', '.conf'
}

# --- РЕГУЛЯРНЫЕ ВЫРАЖЕНИЯ ПАТЧЕРА ---
COMMIT_PATTERN = re.compile(r'<{7}\s*COMMIT\s*\n(.*?)\n>{7}\s*COMMIT', re.DOTALL)
CREATE_PATTERN = re.compile(r'<{7}\s*CREATE:\s*([^\n\r]+?)\s*\n(.*?)\n>{7}\s*CREATE', re.DOTALL)
DELETE_PATTERN = re.compile(r'<{7}\s*DELETE:\s*([^\n\r]+?)\s*\n>{7}\s*DELETE', re.DOTALL)
REPLACE_PATTERN = re.compile(
    r'(?:^|\n)([^\n\r`<>]+?)\n'
    r'<{7}\s*SEARCH\s*\n'
    r'(.*?)\n'
    r'={7}\s*\n'
    r'(.*?)\n'
    r'>{7}\s*REPLACE',
    re.DOTALL
)

# ============================================================================
#                               GIT УТИЛИТЫ
# ============================================================================

def run_git(args, check=True):
    return subprocess.run(["git"] + args, capture_output=True, text=True, check=check)

def check_git_clean():
    """Проверяет, что репозиторий инициализирован и нет незакоммиченных изменений."""
    try:
        res = run_git(["status", "--porcelain"])
        if res.stdout.strip():
            print("❌ Ошибка: В рабочей директории есть незакоммиченные изменения!")
            print("Закоммитьте или спрячьте их (`git stash`), чтобы в случае чего можно было чисто откатиться.")
            sys.exit(1)
    except FileNotFoundError:
        print("❌ Ошибка: Git не установлен или недоступен в PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError:
        print("❌ Ошибка: Текущая директория не является Git-репозиторием. Выполните `git init`.")
        sys.exit(1)

def git_commit(commit_msg):
    run_git(["add", "-A"])
    full_msg = f"[LLM] {commit_msg}"
    run_git(["commit", "-m", full_msg])
    print(f"📦 Сделан Git-коммит: \"{full_msg}\"")

def git_rollback():
    """Откатывает последний коммит, сделанный нейросетью."""
    try:
        # Получаем сообщение последнего коммита
        last_commit = run_git(["log", "-1", "--pretty=%B"]).stdout.strip()
        if not last_commit.startswith("[LLM]"):
            confirm = input(f"⚠️ Последний коммит \"{last_commit}\" не помечен префиксом [LLM]. Точно откатить? (y/N): ")
            if confirm.lower() != 'y':
                print("Откат отменен.")
                return
        
        # Сброс коммита и удаление неотслеживаемых файлов
        run_git(["reset", "--hard", "HEAD~1"])
        run_git(["clean", "-fd"])
        print(f"✨ Успешный откат! Коммит \"{last_commit}\" удален, рабочая копия сброшена.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Ошибка Git при откате: {e.stderr}")

# ============================================================================
#                          ФОРМАТИРОВАНИЕ (BLACK)
# ============================================================================

def format_with_black(target_files):
    py_files = [f for f in target_files if f.endswith(".py") and os.path.exists(f)]
    if not py_files:
        return

    print("🎨 Форматирование кода через black...")
    try:
        subprocess.run([sys.executable, "-m", "black", "-q"] + py_files, check=True)
        print("✅ Код отформатирован.")
    except subprocess.CalledProcessError:
        print("⚠️ Black завершился с предупреждением (возможно синтаксическая ошибка в коде от LLM).")
    except FileNotFoundError:
        print("⚠️ Пакет 'black' не найден. Установите через `pip install black` для автоформатирования.")

# ============================================================================
#                             ПАРСЕР И ПАТЧЕР
# ============================================================================

def apply_patch(patch_text):
    check_git_clean()
    patch_text = patch_text.replace("\r\n", "\n")

    # 1. Извлекаем коммит-сообщение
    commit_match = COMMIT_PATTERN.search(patch_text)
    commit_msg = commit_match.group(1).strip() if commit_match else "Apply LLM changes"

    creates = CREATE_PATTERN.findall(patch_text)
    deletes = DELETE_PATTERN.findall(patch_text)
    replaces = REPLACE_PATTERN.findall(patch_text)

    if not (creates or deletes or replaces):
        print("❌ В переданном тексте не найдено команд (CREATE, DELETE или SEARCH/REPLACE).")
        return

    touched_files = set()

    # 2. Создание файлов
    for file_path, content in creates:
        file_path = file_path.strip()
        path_obj = Path(file_path)
        path_obj.parent.mkdir(parents=True, exist_ok=True)
        with open(path_obj, "w", encoding="utf-8") as f:
            f.write(content.strip() + "\n")
        touched_files.add(file_path)
        print(f"📄 Создан файл: {file_path}")

    # 3. Удаление файлов
    for file_path in deletes:
        file_path = file_path.strip()
        path_obj = Path(file_path)
        if path_obj.exists():
            path_obj.unlink()
            print(f"🗑️ Удален файл: {file_path}")
        else:
            print(f"⚠️ Файл для удаления не найден: {file_path}")

    # 4. Модификация (SEARCH / REPLACE)
    files_data = {}
    replace_targets = set(f.strip() for f, _, _ in replaces)
    
    # Предзагрузка файлов
    for f in replace_targets:
        if not os.path.exists(f):
            print(f"❌ Файл не существует для замены: {f}")
            print("🛑 Откат к исходному состоянию...")
            run_git(["restore", "."])
            run_git(["clean", "-fd"])
            return
        with open(f, "r", encoding="utf-8") as file_obj:
            files_data[f] = file_obj.read()

    # Применение замен
    for file_path, search_block, replace_block in replaces:
        file_path = file_path.strip()
        content = files_data[file_path]

        if search_block not in content:
            print(f"❌ Не найден блок поиска в {file_path}:\n---\n{search_block}\n---")
            run_git(["restore", "."])
            run_git(["clean", "-fd"])
            return

        if content.count(search_block) > 1:
            print(f"❌ Блок поиска в {file_path} не уникален (встречается несколько раз).")
            run_git(["restore", "."])
            run_git(["clean", "-fd"])
            return

        files_data[file_path] = content.replace(search_block, replace_block, 1)
        touched_files.add(file_path)

    # Запись измененных файлов
    for file_path, content in files_data.items():
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"✏️ Изменен файл: {file_path}")

    # 5. Форматирование black
    format_with_black(list(touched_files))

    # 6. Фиксация в Git
    git_commit(commit_msg)

# ============================================================================
#                            ДАМП ПРОЕКТА
# ============================================================================

def is_text_file(file_path: Path) -> bool:
    if file_path.suffix.lower() in TEXT_EXTENSIONS or file_path.name.startswith('.'):
        return True
    try:
        with open(file_path, 'rb') as f:
            chunk = f.read(1024)
            if b'\x00' in chunk:
                return False
            chunk.decode('utf-8')
            return True
    except (UnicodeDecodeError, PermissionError, OSError):
        return False

def build_tree(dir_path: Path, ignore_dirs: set, ignore_files: set, prefix: str = "") -> list:
    tree_lines = []
    try:
        entries = sorted(list(dir_path.iterdir()), key=lambda x: (x.is_file(), x.name.lower()))
    except (PermissionError, OSError):
        return [prefix + "└── [Доступ запрещен]"]

    entries = [
        e for e in entries 
        if (e.is_dir() and e.name not in ignore_dirs) or 
           (e.is_file() and e.name not in ignore_files)
    ]

    total = len(entries)
    for index, entry in enumerate(entries):
        is_last = (index == total - 1)
        connector = "└── " if is_last else "├── "
        
        if entry.is_dir():
            tree_lines.append(f"{prefix}{connector}{entry.name}/")
            extension_prefix = "    " if is_last else "│   "
            tree_lines.extend(build_tree(entry, ignore_dirs, ignore_files, prefix + extension_prefix))
        else:
            tree_lines.append(f"{prefix}{connector}{entry.name}")
            
    return tree_lines

def collect_files(dir_path: Path, ignore_dirs: set, ignore_files: set) -> list:
    files_list = []
    for root, dirs, files in os.walk(dir_path):
        dirs[:] = [d for d in dirs if d not in ignore_dirs]
        for file in files:
            if file in ignore_files:
                continue
            files_list.append(Path(root) / file)
    return sorted(files_list)

def dump_project(root_dir: str, output_file: str):
    root_path = Path(root_dir).resolve()
    out_path = Path(output_file).resolve()
    
    ignore_dirs = DEFAULT_IGNORE_DIRS.copy()
    ignore_files = DEFAULT_IGNORE_FILES.copy()
    ignore_files.add(out_path.name)

    try:
        with open(out_path, 'w', encoding='utf-8') as out:
            out.write("================================================\n")
            out.write(f" СТРУКТУРА ПРОЕКТА: {root_path.name}\n")
            out.write("================================================\n\n")
            out.write(f"{root_path.name}/\n")
            
            tree_lines = build_tree(root_path, ignore_dirs, ignore_files)
            out.write("\n".join(tree_lines))
            out.write("\n\n\n")

            out.write("================================================\n")
            out.write(" СОДЕРЖИМОЕ ФАЙЛОВ\n")
            out.write("================================================\n\n")

            all_files = collect_files(root_path, ignore_dirs, ignore_files)
            
            for file_path in all_files:
                try:
                    rel_path = file_path.relative_to(root_path)
                except ValueError:
                    rel_path = file_path
                
                if not is_text_file(file_path):
                    continue
                    
                out.write(f"\n{'='*70}\n")
                out.write(f"FILE: {rel_path}\n")
                out.write(f"{'='*70}\n")
                
                try:
                    with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
                        content = f.read()
                        out.write(content)
                        if not content.endswith('\n'):
                            out.write('\n')
                except (PermissionError, OSError) as e:
                    out.write(f"[Ошибка доступа к файлу: {e}]\n")

        print(f"✅ Дамп готов! Результат сохранен в: {out_path}")
    except Exception as e:
        print(f"❌ Критическая ошибка записи дампа: {e}")

# ============================================================================
#                                CLI ИНТЕРФЕЙС
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="AI DevTool: Дамп, автопатчинг через Git и black.")
    subparsers = parser.add_subparsers(dest="command")

    # Команда дампа
    dump_parser = subparsers.add_parser("dump", help="Собрать дамп проекта в один файл")
    dump_parser.add_argument("path", nargs="?", default=".", help="Путь к проекту (по умолчанию .)")
    dump_parser.add_argument("-o", "--output", default="project_dump.txt", help="Выходной файл")

    # Команда применения
    apply_parser = subparsers.add_parser("apply", help="Применить правки от нейросети")
    apply_parser.add_argument("file", nargs="?", help="Файл с ответом нейросети (если не указан, читает stdin)")

    # Команда отката
    subparsers.add_parser("rollback", help="Откатить последний коммит правок [LLM]")

    args = parser.parse_args()

    if args.command == "dump":
        dump_project(args.path, args.output)

    elif args.command == "apply":
        if args.file:
            with open(args.file, "r", encoding="utf-8") as f:
                patch_text = f.read()
        else:
            print("Вставьте ответ нейросети и нажмите Ctrl+D (Linux/macOS) или Ctrl+Z + Enter (Windows):")
            patch_text = sys.stdin.read()
        apply_patch(patch_text)

    elif args.command == "rollback":
        git_rollback()

    else:
        parser.print_help()

if __name__ == "__main__":
    main()