param(
    [Parameter(Mandatory = $true)]
    [string[]] $PredComplex,
    [string] $Source = "boltz"
)

$ErrorActionPreference = "Stop"
$ModeRoot = Split-Path -Parent $PSScriptRoot
$PipelineRoot = Split-Path -Parent $ModeRoot
$PyMol = "C:\ProgramData\pymol\Scripts\pymol.exe"
$AlignScript = Join-Path $PipelineRoot "05b_align_pose.pml"
$OutputDir = Join-Path $ModeRoot "cofold_poses"

if (-not (Test-Path -LiteralPath $PyMol)) {
    throw "PyMOL executable not found: $PyMol"
}
if (-not (Test-Path -LiteralPath $AlignScript)) {
    throw "Alignment script not found: $AlignScript"
}

$templates = @(
    @{ Name = "9WZU_apo"; Path = Join-Path $PipelineRoot "templates\9WZU_apo.pdb" },
    @{ Name = "apo_8WL6"; Path = Join-Path $PipelineRoot "templates\apo_8WL6.pdb" },
    @{ Name = "caspo_T2"; Path = Join-Path $PipelineRoot "templates\caspo_T2.pdb" }
)

New-Item -ItemType Directory -Path $OutputDir -Force | Out-Null
foreach ($complex in $PredComplex) {
    $resolvedComplex = (Resolve-Path -LiteralPath $complex).Path
    $sample = [IO.Path]::GetFileNameWithoutExtension($resolvedComplex) -replace '[^A-Za-z0-9_.-]', '_'
    foreach ($template in $templates) {
        if (-not (Test-Path -LiteralPath $template.Path)) {
            throw "Template not found: $($template.Path)"
        }
        $out = Join-Path $OutputDir ("{0}__{1}_{2}.sdf" -f $template.Name, $Source, $sample)
        Write-Host "[align] $sample -> $($template.Name)"
        & $PyMol -cq $AlignScript -- $template.Path $resolvedComplex $out
        if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $out)) {
            throw "PyMOL alignment failed for $sample -> $($template.Name)"
        }
    }
}

Write-Host "Aligned poses written to: $OutputDir"
