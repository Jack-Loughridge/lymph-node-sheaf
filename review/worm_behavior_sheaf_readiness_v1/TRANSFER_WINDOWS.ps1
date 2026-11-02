# No full-collection download or raw microscopy. Download the five named files
# in your browser first; then run this script in Windows PowerShell.
$InputFolder = Read-Host 'Full Windows folder containing the five Dryad files'
$VmDir = '/home/Jack/Real_Data/lymph_node_sheaf/data/worm_behavior_sheaf_readiness_v1/recordings'
$Files = @('2023-03-30-01-cudata.h5','2023-06-21-01-cudata.h5','2023-06-23-08-cudata.h5','NeuroPAL_labels.jld2','load_data_example.jl')
$Receipt = @()
ssh tda-vm "mkdir -p '$VmDir'"
if ($LASTEXITCODE -ne 0) { throw 'Cannot prepare VM destination' }
foreach ($Name in $Files) {
    $Source = Join-Path $InputFolder $Name
    if (-not (Test-Path -LiteralPath $Source -PathType Leaf)) { throw "Missing $Name" }
    $Hash = (Get-FileHash -LiteralPath $Source -Algorithm SHA256).Hash.ToLowerInvariant()
    $Size = (Get-Item -LiteralPath $Source).Length
    ssh tda-vm "test ! -e '$VmDir/$Name'"
    if ($LASTEXITCODE -ne 0) { throw "Destination exists; inspect instead of overwrite: $Name" }
    scp -- "$Source" "tda-vm:$VmDir/"
    if ($LASTEXITCODE -ne 0) { throw "SCP failed: $Name" }
    $Remote = ssh tda-vm "sha256sum '$VmDir/$Name'"
    if ($LASTEXITCODE -ne 0) { throw "Remote checksum failed: $Name" }
    $RemoteHash = ($Remote -split '\s+')[0].ToLowerInvariant()
    if ($RemoteHash -ne $Hash) { throw "Checksum mismatch: $Name; preserve for inspection" }
    $Receipt += [pscustomobject]@{file=$Name;bytes=$Size;sha256=$Hash;transfer_verified=$true}
    Write-Host "${Name}: checksum verified"
}
$Receipt | ConvertTo-Json | Set-Content -Encoding utf8 (Join-Path $InputFolder 'transfer_receipt_PRIVATE.json')
Write-Host 'Sample transferred; scientific QC is still required. Do not download all 27 yet.'
