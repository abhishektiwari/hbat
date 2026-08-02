"""Fail if changed docs contain invisible/hidden Unicode (zero-width, bidi-override, control chars)."""
import os
import sys
import unicodedata

SUSPICIOUS_CODEPOINTS = {
    0x200B: "ZERO WIDTH SPACE",
    0x200C: "ZERO WIDTH NON-JOINER",
    0x200D: "ZERO WIDTH JOINER",
    0x200E: "LEFT-TO-RIGHT MARK",
    0x200F: "RIGHT-TO-LEFT MARK",
    0x2060: "WORD JOINER",
    0xFEFF: "ZERO WIDTH NO-BREAK SPACE (BOM)",
    0x202A: "LEFT-TO-RIGHT EMBEDDING",
    0x202B: "RIGHT-TO-LEFT EMBEDDING",
    0x202C: "POP DIRECTIONAL FORMATTING",
    0x202D: "LEFT-TO-RIGHT OVERRIDE",
    0x202E: "RIGHT-TO-LEFT OVERRIDE",
    0x2066: "LEFT-TO-RIGHT ISOLATE",
    0x2067: "RIGHT-TO-LEFT ISOLATE",
    0x2068: "FIRST STRONG ISOLATE",
    0x2069: "POP DIRECTIONAL ISOLATE",
}


def scan(path):
    try:
        with open(path, encoding="utf-8") as f:
            content = f.read()
    except (FileNotFoundError, UnicodeDecodeError):
        return []

    findings = []
    for i, ch in enumerate(content):
        cp = ord(ch)
        is_suspicious = cp in SUSPICIOUS_CODEPOINTS
        is_control = (0x0000 <= cp <= 0x001F and ch not in "\n\t\r") or (0x007F <= cp <= 0x009F)
        if is_suspicious or is_control:
            name = SUSPICIOUS_CODEPOINTS.get(cp, unicodedata.name(ch, "UNKNOWN CONTROL CHAR"))
            line_no = content.count("\n", 0, i) + 1
            findings.append(f"- `{path}` line {line_no}: U+{cp:04X} ({name})")
    return findings


def main():
    files = sys.argv[1:]
    findings = []
    for path in files:
        findings.extend(scan(path))

    if not findings:
        print("No hidden/invisible Unicode characters found.")
        return 0

    print("::error::Hidden/invisible Unicode characters found in changed documentation — possible prompt-injection or homoglyph attack")

    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary_path:
        with open(summary_path, "a") as f:
            f.write("## Hidden Unicode characters found in documentation\n\n")
            f.write(
                "These characters are invisible when rendered but can hide instructions from "
                "AI tools or humans (e.g. prompt injection, bidi-override homoglyph attacks). "
                "Manual review required.\n\n"
            )
            f.write("\n".join(findings) + "\n")

    return 1


if __name__ == "__main__":
    sys.exit(main())
