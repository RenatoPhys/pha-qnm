$ErrorActionPreference = 'Stop'

$paperRoot = $PSScriptRoot
$toolRoot = Join-Path $paperRoot 'tools\tectonic'
$tectonic = Join-Path $toolRoot 'tectonic.exe'
$archive = Join-Path $paperRoot 'tools\tectonic.zip'
$version = '0.17.0'
$expectedSha256 = 'F61CE51F0B0ADE1015B7DE7EF368541C5424E9756ECBD0D7AF97D6D48030845F'
$url = "https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%40$version/tectonic-$version-x86_64-pc-windows-msvc.zip"

if (-not (Test-Path -LiteralPath $tectonic)) {
    New-Item -ItemType Directory -Force -Path $toolRoot | Out-Null
    Invoke-WebRequest -Uri $url -OutFile $archive
    $actualSha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $archive).Hash
    if ($actualSha256 -ne $expectedSha256) {
        Remove-Item -LiteralPath $archive -Force
        throw "Tectonic archive checksum mismatch: $actualSha256"
    }
    Expand-Archive -LiteralPath $archive -DestinationPath $toolRoot -Force
    Remove-Item -LiteralPath $archive -Force
}

Push-Location $paperRoot
try {
    & $tectonic '.\main.tex' --keep-logs --keep-intermediates
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} finally {
    Pop-Location
}

Write-Output "PDF generated: $(Join-Path $paperRoot 'main.pdf')"
