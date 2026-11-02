# Transfer only the remaining 24 public processed recordings. No raw microscopy.
$InputFolder = Read-Host 'Windows directory containing the downloaded/extracted remaining HDF5 files'
$VmDir = '/home/Jack/Real_Data/lymph_node_sheaf/data/worm_behavior_sheaf_readiness_v1/recordings'
$VmPrivate = '/home/Jack/Real_Data/lymph_node_sheaf/data/worm_behavior_sheaf_readiness_v1/export_recovery_v3'
$Names = @(
    '2023-06-29-01-cudata.h5',
    '2023-06-29-13-cudata.h5',
    '2023-07-14-08-cudata.h5',
    '2023-07-14-14-cudata.h5',
    '2023-07-27-01-cudata.h5',
    '2023-07-27-08-cudata.h5',
    '2023-08-08-07-cudata.h5',
    '2023-08-14-01-cudata.h5',
    '2023-08-15-01-cudata.h5',
    '2023-08-16-01-cudata.h5',
    '2023-08-21-01-cudata.h5',
    '2023-09-07-01-cudata.h5',
    '2023-09-14-01-cudata.h5',
    '2023-09-25-01-cudata.h5',
    '2023-10-05-01-cudata.h5',
    '2023-12-11-01-cudata.h5',
    '2024-02-01-07-cudata.h5',
    '2024-05-22-15-cudata.h5',
    '2024-05-22-21-cudata.h5',
    '2024-05-28-02-cudata.h5',
    '2024-05-28-10-cudata.h5',
    '2024-06-18-01-cudata.h5',
    '2024-06-18-07-cudata.h5',
    '2024-06-18-12-cudata.h5'
)
$Stamp = [DateTime]::UtcNow.ToString('yyyyMMddTHHmmssZ')
$ReceiptName = "remaining_transfer_${Stamp}_PRIVATE.json"
$ReceiptPath = Join-Path $InputFolder $ReceiptName
$Receipt = @()
$Sources = Get-ChildItem -LiteralPath $InputFolder -Recurse -File
ssh tda-vm "mkdir -p '$VmDir' '$VmPrivate'"
if ($LASTEXITCODE -ne 0) { throw 'VM directory setup failed' }
foreach ($Name in $Names) {
    $Match = @($Sources | Where-Object { $_.Name -eq $Name })
    if ($Match.Count -ne 1) { throw "Expected exactly one local file named $Name; found $($Match.Count)" }
    $Source = $Match[0]
    if ($Source.Length -gt 20000000) { throw "Unexpected processed file size: $Name" }
    $Hash = (Get-FileHash -LiteralPath $Source.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    ssh tda-vm "test ! -e '$VmDir/$Name'"
    if ($LASTEXITCODE -ne 0) { throw "Destination exists; preserve and inspect: $Name" }
    scp -- "$($Source.FullName)" "tda-vm:$VmDir/"
    if ($LASTEXITCODE -ne 0) { throw "SCP failed: $Name" }
    $Remote = ssh tda-vm "sha256sum '$VmDir/$Name'"
    if ($LASTEXITCODE -ne 0) { throw "Remote hash failed: $Name" }
    if (($Remote -split '\s+')[0].ToLowerInvariant() -ne $Hash) { throw "Transfer hash mismatch: $Name" }
    $Receipt += [pscustomobject]@{filename=$Name;bytes=$Source.Length;sha256=$Hash}
    @{algorithm='SHA256';files=$Receipt} | ConvertTo-Json -Depth 4 | Set-Content -Encoding utf8 $ReceiptPath
    Write-Host "${Name}: SHA-256 verified; Windows receipt checkpointed"
}
scp -- "$ReceiptPath" "tda-vm:$VmPrivate/"
if ($LASTEXITCODE -ne 0) { throw 'Recordings transferred but private receipt copy failed; preserve the Windows receipt' }
Write-Host "Complete. Private VM receipt: $VmPrivate/$ReceiptName"
