# JHEP submission checklist

Checked against the official JHEP author instructions and tool downloads on
2026-08-25.

## Closed locally

- `paper/jheppub.sty` is byte-identical to the official JHEP download
  (SHA-256 `478C47BF3D831B723766104B4C506DE19F04BF5EDE537CC25AEDA9F3F5ADE2A0`).
- `paper/JHEP.bst` is byte-identical to the official version 2.18
  (SHA-256 `7D94537265BA1185F58D85B052D92BAC855A0740DDDDBCDBE338621CC4649F4A`).
- The abstract fits on page 1 and contains neither citations nor displayed
  formulae.
- The four manuscript keywords are exact entries from the controlled JHEP
  keyword list.
- The acknowledgments disclose the uses of OpenAI Codex in preparing the
  manuscript, as required by the current AI-assisted-technology policy.
- The source build generates `main.bbl`; the submission packager requires it.
- The master file is `main.tex` at the archive root, with figures under
  `figures/`.
- Only compilation inputs enter the source archive: `main.tex`, `main.bbl`,
  `references.bib`, `jheppub.sty`, `JHEP.bst`, and the eight referenced PDF
  figures.
- The 13-file archive layout was assembled and compiled successfully in an
  isolated directory with no repository files available implicitly.
- All figure and manuscript PDF fonts are embedded, and no Type 3 fonts remain.
  The manuscript compiles
  without undefined citations/references, duplicate labels, overfull boxes, or
  TeX errors.
- The reproducibility appendix identifies the intended code repository and the
  upstream posterior-data citation without claiming that the still-local
  reconstruction has already been released.

## Author-controlled items still required

- Insert every author's first and family names, affiliation, and email in
  `main.tex`; designate the corresponding author.
- Replace the author marker in `CITATION.cff` and supply ORCID identifiers.
- Add funding identifiers and any additional acknowledgments, or explicitly
  confirm that there are none.
- Obtain the arXiv identifier, add `\arxivnumber{...}`, and ensure that the
  submitted JHEP version is identical to the arXiv version.
- Approve authorship, originality, institutional authorization, and exclusive
  submission declarations.
- Select two to four of the manuscript's controlled keywords in the JHEP form
  and choose the matching data/code-availability statement.
- Publish the reconstructed branch and machine-readable outputs at the named
  repository no later than submission.
- Archive the code release in a persistent repository and provide its DOI when
  requested at proofreading.

## Final packaging command

```powershell
powershell -ExecutionPolicy Bypass -File paper/build_submission_archive.ps1
```

The command deliberately fails while author, affiliation, email, arXiv, or
`CITATION.cff` authorship metadata are missing. On success it emits only the
necessary compilation files as `paper/pha-qnm-jhep-source.tar.gz` and prints
the archive SHA-256 digest.

Official sources:

- <https://jhep.sissa.it/jhep/help/JHEP/JHEP_author.jsp>
- <https://jhep.sissa.it/jhep/help/JHEP_TeXclass.jsp>
- <https://jhep.sissa.it/jhep/help/keywordsList.jsp>
