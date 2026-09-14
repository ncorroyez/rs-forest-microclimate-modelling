#!/usr/bin/env python3
"""Build the Annals of Forest Science submission reference.docx for pandoc.

Patches pandoc's own default reference.docx into the journal template:
Cambria 12 pt body, Calibri headings, double spacing, continuous line
numbering, A4 with 2.5 cm margins, and a "Page X of Y" footer.

    python3 Chapitre1/redaction/tools/build_aofs_reference.py [outdir]

Writes <outdir>/reference_aofs.docx (default: alongside this script).
It lived in /tmp twice and evaporated twice; keep it in the repository.

Every substitution is asserted. A patch that silently matches nothing used
to produce an unpatched file that still looked like a success, so rw()
raises instead.
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE
BUILD = os.path.join(OUTDIR, "_aofs_build")
X = os.path.join(BUILD, "x")

BODY_FONT = "Cambria"
HEAD_FONT = "Calibri"
# 480 twips = 24 pt = double the 12 pt body
LINE = '<w:spacing w:after="200" w:line="480" w:lineRule="auto" />'
# A4 in twips; 1418 twips = 2.5 cm
SECTPR = (
    "<w:sectPr>"
    '<w:footerReference w:type="default" r:id="rIdFooterAoFS" />'
    '<w:pgSz w:w="11906" w:h="16838" />'
    '<w:pgMar w:top="1418" w:right="1418" w:bottom="1418" w:left="1418"'
    ' w:header="708" w:footer="708" w:gutter="0" />'
    '<w:lnNumType w:countBy="1" w:restart="continuous" w:distance="360" />'
    "</w:sectPr>"
)
FOOTER = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
    '<w:p><w:pPr><w:jc w:val="center" /></w:pPr>'
    "<w:r><w:t xml:space=\"preserve\">Page </w:t></w:r>"
    '<w:r><w:fldChar w:fldCharType="begin" /></w:r>'
    '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
    '<w:r><w:fldChar w:fldCharType="separate" /></w:r>'
    "<w:r><w:t>1</w:t></w:r>"
    '<w:r><w:fldChar w:fldCharType="end" /></w:r>'
    "<w:r><w:t xml:space=\"preserve\"> of </w:t></w:r>"
    '<w:r><w:fldChar w:fldCharType="begin" /></w:r>'
    '<w:r><w:instrText xml:space="preserve"> NUMPAGES </w:instrText></w:r>'
    '<w:r><w:fldChar w:fldCharType="separate" /></w:r>'
    "<w:r><w:t>1</w:t></w:r>"
    '<w:r><w:fldChar w:fldCharType="end" /></w:r>'
    "</w:p></w:ftr>"
)


def rw(path, fn):
    """Rewrite a file inside the extracted docx, failing loudly on a no-op."""
    p = os.path.join(X, path)
    with open(p, encoding="utf8") as fh:
        s = fh.read()
    s2 = fn(s)
    assert s2 != s, f"no change made to {path}"
    with open(p, "w", encoding="utf8") as fh:
        fh.write(s2)


def main():
    shutil.rmtree(BUILD, ignore_errors=True)
    os.makedirs(X)
    base = os.path.join(BUILD, "reference.docx")
    with open(base, "wb") as fh:
        fh.write(subprocess.run(
            ["pandoc", "--print-default-data-file", "reference.docx"],
            check=True, stdout=subprocess.PIPE).stdout)
    subprocess.run(["unzip", "-q", base, "-d", X], check=True)

    # 1. fonts and double spacing in docDefaults, Calibri on the heading styles
    def styles(s):
        s = re.sub(
            r'<w:rFonts w:asciiTheme="minorHAnsi"[^/]*/>',
            f'<w:rFonts w:ascii="{BODY_FONT}" w:eastAsia="{BODY_FONT}"'
            f' w:hAnsi="{BODY_FONT}" w:cs="{BODY_FONT}" />', s, count=1)
        s = s.replace('<w:spacing w:after="200" />', LINE, 1)

        def head(m):
            blk = m.group(0)
            if not re.search(r'w:styleId="(Heading[1-6]|Title|Subtitle)"', blk):
                return blk
            if "<w:rPr>" not in blk:
                return blk
            return blk.replace(
                "<w:rPr>",
                f'<w:rPr><w:rFonts w:ascii="{HEAD_FONT}" w:eastAsia="{HEAD_FONT}"'
                f' w:hAnsi="{HEAD_FONT}" w:cs="{HEAD_FONT}" />', 1)
        return re.sub(r"<w:style [^>]*>.*?</w:style>", head, s, flags=re.S)

    rw("word/styles.xml", styles)

    # 2. page setup, line numbering, footer reference
    rw("word/document.xml",
       lambda s: re.sub(r"<w:sectPr>.*?</w:sectPr>", SECTPR, s, count=1, flags=re.S))

    # 3. the footer part itself, its relationship and its content type
    with open(os.path.join(X, "word/footer1.xml"), "w", encoding="utf8") as fh:
        fh.write(FOOTER)
    rw("word/_rels/document.xml.rels", lambda s: s.replace(
        "</Relationships>",
        '<Relationship Id="rIdFooterAoFS"'
        ' Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer"'
        ' Target="footer1.xml" /></Relationships>', 1))
    rw("[Content_Types].xml", lambda s: s.replace(
        "</Types>",
        '<Override PartName="/word/footer1.xml"'
        ' ContentType="application/vnd.openxmlformats-officedocument'
        '.wordprocessingml.footer+xml" /></Types>', 1))

    out = os.path.join(OUTDIR, "reference_aofs.docx")
    if os.path.exists(out):
        os.remove(out)
    subprocess.run(["zip", "-q", "-r", "-X", out, "."], cwd=X, check=True)
    shutil.rmtree(BUILD, ignore_errors=True)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
