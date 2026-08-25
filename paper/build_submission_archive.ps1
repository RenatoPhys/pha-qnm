param(
    [string]$Output = 'pha-qnm-jhep-source.tar.gz'
)

$ErrorActionPreference = 'Stop'

$paperRoot = [System.IO.Path]::GetFullPath($PSScriptRoot)
$repoRoot = [System.IO.Path]::GetFullPath((Join-Path $paperRoot '..'))
$mainTex = Join-Path $paperRoot 'main.tex'
$citationFile = Join-Path $repoRoot 'CITATION.cff'
$source = Get-Content -LiteralPath $mainTex -Raw
$citation = Get-Content -LiteralPath $citationFile -Raw

$requiredMetadata = [ordered]@{
    author = '\\author\s*\{'
    affiliation = '\\affiliation(?:\[[^]]*\])?\s*\{'
    email = '\\emailAdd\s*\{'
    arXiv = '\\arxivnumber\s*\{'
}

$missingMetadata = @()
foreach ($entry in $requiredMetadata.GetEnumerator()) {
    if ($source -notmatch $entry.Value) {
        $missingMetadata += $entry.Key
    }
}
if ($citation -match '\[AUTHOR TO BE SUPPLIED\]') {
    $missingMetadata += 'CITATION.cff author'
}
if ($missingMetadata.Count -gt 0) {
    throw ('Submission metadata are incomplete: ' + ($missingMetadata -join ', '))
}

$requiredFiles = @(
    'main.tex',
    'main.bbl',
    'references.bib',
    'jheppub.sty',
    'JHEP.bst'
)

$figureMatches = [regex]::Matches(
    $source,
    '\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}'
)
$figureFiles = @($figureMatches | ForEach-Object { $_.Groups[1].Value } | Sort-Object -Unique)
$archiveFiles = @($requiredFiles + $figureFiles)

foreach ($relativePath in $archiveFiles) {
    $absolutePath = Join-Path $paperRoot $relativePath
    if (-not (Test-Path -LiteralPath $absolutePath -PathType Leaf)) {
        throw "Required submission file is missing: $relativePath"
    }
}

$paperTmp = [System.IO.Path]::GetFullPath((Join-Path $paperRoot 'tmp'))
New-Item -ItemType Directory -Force -Path $paperTmp | Out-Null
$stage = [System.IO.Path]::GetFullPath((Join-Path $paperTmp ('submission-' + [guid]::NewGuid().ToString('N'))))
if (-not $stage.StartsWith($paperTmp + [System.IO.Path]::DirectorySeparatorChar, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Unsafe staging path: $stage"
}
New-Item -ItemType Directory -Path $stage | Out-Null

try {
    foreach ($relativePath in $archiveFiles) {
        $destination = Join-Path $stage $relativePath
        $destinationDirectory = Split-Path -Parent $destination
        if ($destinationDirectory) {
            New-Item -ItemType Directory -Force -Path $destinationDirectory | Out-Null
        }
        Copy-Item -LiteralPath (Join-Path $paperRoot $relativePath) -Destination $destination
    }

    $archive = [System.IO.Path]::GetFullPath((Join-Path $paperRoot $Output))
    if (-not $archive.StartsWith($paperRoot + [System.IO.Path]::DirectorySeparatorChar, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Archive output must remain under the paper directory: $archive"
    }
    if (Test-Path -LiteralPath $archive) {
        Remove-Item -LiteralPath $archive -Force
    }

    & tar.exe -czf $archive -C $stage @archiveFiles
    if ($LASTEXITCODE -ne 0) {
        throw 'tar failed while creating the JHEP source archive'
    }

    $listed = @(& tar.exe -tzf $archive)
    if ($LASTEXITCODE -ne 0) {
        throw 'tar failed while validating the JHEP source archive'
    }
    $unexpected = @($listed | Where-Object { $_ -notin $archiveFiles })
    $missing = @($archiveFiles | Where-Object { $_ -notin $listed })
    if ($unexpected.Count -gt 0 -or $missing.Count -gt 0) {
        throw "Archive content mismatch. Missing: $($missing -join ', '); unexpected: $($unexpected -join ', ')"
    }

    $hash = (Get-FileHash -Algorithm SHA256 -LiteralPath $archive).Hash
    Write-Output "JHEP source archive: $archive"
    Write-Output "Files: $($archiveFiles.Count)"
    Write-Output "SHA256: $hash"
}
finally {
    if (Test-Path -LiteralPath $stage) {
        Remove-Item -LiteralPath $stage -Recurse -Force
    }
}
