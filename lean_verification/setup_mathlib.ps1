# Robust Mathlib setup for Windows
# Run from PowerShell: .\setup_mathlib.ps1
# Or: powershell -ExecutionPolicy Bypass -File setup_mathlib.ps1

$env:PATH = "$env:USERPROFILE\.elan\bin;$env:PATH"

Write-Host "============================================"
Write-Host "Mathlib Setup for TwistorVerified"
Write-Host "============================================"
Write-Host ""

Set-Location $PSScriptRoot

# Step 1: lake update (clones Mathlib)
$maxRetries = 10
$retryDelay = 30

for ($i = 1; $i -le $maxRetries; $i++) {
    Write-Host "[Attempt $i/$maxRetries] Running lake update..."
    $result = & lake update 2>&1
    $exitCode = $LASTEXITCODE
    Write-Host $result

    if ($exitCode -eq 0) {
        Write-Host "lake update succeeded!"
        break
    }

    if ($i -lt $maxRetries) {
        Write-Host "Failed (exit $exitCode). Retrying in ${retryDelay}s..."
        Start-Sleep -Seconds $retryDelay
    }
}

if ($exitCode -ne 0) {
    Write-Host "lake update failed after $maxRetries attempts."
    exit 1
}

# Step 2: Download prebuilt oleans
Write-Host ""
Write-Host "Downloading prebuilt Mathlib cache (~4GB)..."
Write-Host "This avoids building Mathlib from source (which takes hours)."

for ($i = 1; $i -le $maxRetries; $i++) {
    Write-Host "[Attempt $i/$maxRetries] Running lake exe cache get..."
    $result = & lake exe cache get 2>&1
    $exitCode = $LASTEXITCODE
    Write-Host $result

    if ($exitCode -eq 0) {
        Write-Host "Cache download succeeded!"
        break
    }

    if ($i -lt $maxRetries) {
        Write-Host "Failed (exit $exitCode). Retrying in ${retryDelay}s..."
        Start-Sleep -Seconds $retryDelay
    }
}

# Step 3: Build
Write-Host ""
Write-Host "Building TwistorVerified..."
& lake build 2>&1

if ($LASTEXITCODE -eq 0) {
    Write-Host ""
    Write-Host "============================================"
    Write-Host "SUCCESS! Mathlib is ready."
    Write-Host "You can now use norm_num, interval_cases, etc."
    Write-Host "============================================"
} else {
    Write-Host ""
    Write-Host "Build had errors - check above."
    Write-Host "The .lean files may need updating for the new Mathlib API."
}
