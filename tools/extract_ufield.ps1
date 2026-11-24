# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later
# Extracts DELPHI's measured UFIELD barrel & endcap lookup tables from ufield.car (UFIDAT)
# into core/src/ufield_data.cpp.
param(
  [string]$UfieldCar = "simana_v94d\src\car\ufield.car",
  [string]$OutCpp    = "delphipp\core\src\ufield_data.cpp"
)

$ErrorActionPreference = 'Stop'
$lines = [IO.File]::ReadAllLines($UfieldCar)
$inDat = $false
$dataMap = @{}
$curName = $null
$buf = New-Object System.Text.StringBuilder

foreach ($line in $lines) {
  if ($line -match '^\+DECK,UFIDAT') { $inDat = $true; continue }
  if (-not $inDat) { continue }
  if ($line.Length -eq 0 -or $line[0] -eq '*' -or $line[0] -eq 'C') { continue }

  if ($line -match '^\s{6}\s*DATA\s+([A-Z0-9_]+)\s*(.*)$') {
    if ($curName) { $dataMap[$curName] = $buf.ToString() }
    $curName = $matches[1]
    [void]$buf.Clear()
    [void]$buf.Append($matches[2])
  } elseif ($curName -and $line.Length -ge 6 -and $line[5] -ne ' ') {
    [void]$buf.Append(' ')
    [void]$buf.Append($line.Substring(6))
  } elseif ($line -match '^\s{6}\s*END\b') {
    if ($curName) { $dataMap[$curName] = $buf.ToString(); $curName = $null }
  }
}
if ($curName) { $dataMap[$curName] = $buf.ToString() }

function Parse-Nums([string]$raw) {
  if (-not $raw) { return @() }
  $m = [regex]::Match($raw, '/([^/]+)/')
  if (-not $m.Success) { return @() }
  $tokens = $m.Groups[1].Value -split ',' | % { $_.Trim() } | ? { $_ -ne '' }
  $out = New-Object System.Collections.Generic.List[string]
  foreach ($t in $tokens) {
    if ($t -match '^(\d+)\*(.+)$') {
      $cnt = [int]$matches[1]; $val = $matches[2].Trim()
      for ($i=0; $i -lt $cnt; $i++) { $out.Add($val) }
    } else {
      $out.Add($t)
    }
  }
  return $out
}

function Concat-Arrays([string[]]$names) {
  $res = New-Object System.Collections.Generic.List[string]
  foreach ($n in $names) {
    foreach ($v in (Parse-Nums $dataMap[$n])) { $res.Add($v) }
  }
  return $res
}

$rgridb = Parse-Nums $dataMap['RGRIDB']
$zgride = Parse-Nums $dataMap['ZGRIDE']
$nrbarr = Parse-Nums $dataMap['NRBARR']
$nzbarr = Parse-Nums $dataMap['NZBARR']
$nrendc = Parse-Nums $dataMap['NRENDC']
$nzendc = Parse-Nums $dataMap['NZENDC']

$zrba = Concat-Arrays @('ZRBA01','ZRBA02','ZRBA03')
$brba = Concat-Arrays @('BRBA01','BRBA02','BRBA03')
$zzba = Concat-Arrays @('ZZBA01','ZZBA02','ZZBA03')
$bzba = Concat-Arrays @('BZBA01','BZBA02','BZBA03')

$rren = Concat-Arrays @('RREN01','RREN02','RREN03','RREN04','RREN05','RREN06')
$bren = Concat-Arrays @('BREN01','BREN02','BREN03','BREN04','BREN05','BREN06')
$rzen = Concat-Arrays @('RZEN01','RZEN02','RZEN03','RZEN04','RZEN05','RZEN06')
$bzen = Concat-Arrays @('BZEN01','BZEN02','BZEN03','BZEN04','BZEN05','BZEN06')

function Fmt-Float([string]$s) {
  $v = [double]::Parse($s, [Globalization.CultureInfo]::InvariantCulture)
  return ($v.ToString("G7", [Globalization.CultureInfo]::InvariantCulture) + "f")
}

$sb = New-Object System.Text.StringBuilder
[void]$sb.AppendLine("// SPDX-FileCopyrightText: 1986-2024 Members of the DELPHI collaboration")
[void]$sb.AppendLine("// SPDX-FileCopyrightText: 2026 DELPHI++ Contributors")
[void]$sb.AppendLine("// SPDX-License-Identifier: GPL-3.0-or-later")
[void]$sb.AppendLine("//")
[void]$sb.AppendLine("// Auto-generated from simana_v94d/src/car/ufield.car (BLOCKDATA UFIDAT).")
[void]$sb.AppendLine("")
[void]$sb.AppendLine('#include "delphipp/magnetic_field.hpp"')
[void]$sb.AppendLine("")
[void]$sb.AppendLine("namespace delphipp {")
[void]$sb.AppendLine("")
[void]$sb.AppendLine("void MagneticField::load_delphi_map() {")
[void]$sb.AppendLine("  set_scale(1.0f);")
[void]$sb.AppendLine("  set_sat_quad_gradient(0.0007609f);")
[void]$sb.AppendLine("  std::vector<FieldLayerTable> barrel;")
[void]$sb.AppendLine("  barrel.reserve(41);")

for ($lay = 0; $lay -lt 41; $lay++) {
  $rmin = Fmt-Float $rgridb[$lay]
  $rmax = Fmt-Float $rgridb[$lay + 1]
  $nz = [int]$nzbarr[$lay]
  $nr = [int]$nrbarr[$lay]
  $zzs = (0..($nz-1) | % { Fmt-Float $zzba[$lay*4 + $_] }) -join ', '
  $bzs = (0..($nz-1) | % { Fmt-Float $bzba[$lay*4 + $_] }) -join ', '
  $zrs = (0..($nr-1) | % { Fmt-Float $zrba[$lay*4 + $_] }) -join ', '
  $brs = (0..($nr-1) | % { Fmt-Float $brba[$lay*4 + $_] }) -join ', '
  [void]$sb.AppendLine("  barrel.push_back({$rmin, $rmax, {$zzs}, {$bzs}, {$zrs}, {$brs}});")
}
[void]$sb.AppendLine("  set_barrel_layers(std::move(barrel));")
[void]$sb.AppendLine("")
[void]$sb.AppendLine("  std::vector<FieldLayerTable> endcap;")
[void]$sb.AppendLine("  endcap.reserve(39);")
# IGEPL1 = 7 (1-based) -> index 6 (0-based) in ZGRIDE(46) -> 39 layers [6..44]..[7..45]
for ($lay = 0; $lay -lt 39; $lay++) {
  $zmin = Fmt-Float $zgride[6 + $lay]
  $zmax = Fmt-Float $zgride[7 + $lay]
  $nz = [int]$nzendc[$lay]
  $nr = [int]$nrendc[$lay]
  $rzs = (0..($nz-1) | % { Fmt-Float $rzen[$lay*10 + $_] }) -join ', '
  $bzs = (0..($nz-1) | % { Fmt-Float $bzen[$lay*10 + $_] }) -join ', '
  $rrs = (0..($nr-1) | % { Fmt-Float $rren[$lay*10 + $_] }) -join ', '
  $brs = (0..($nr-1) | % { Fmt-Float $bren[$lay*10 + $_] }) -join ', '
  [void]$sb.AppendLine("  endcap.push_back({$zmin, $zmax, {$rzs}, {$bzs}, {$rrs}, {$brs}});")
}
[void]$sb.AppendLine("  set_endcap_layers(std::move(endcap));")
[void]$sb.AppendLine("}")
[void]$sb.AppendLine("")
[void]$sb.AppendLine("}  // namespace delphipp")

[IO.File]::WriteAllText($OutCpp, $sb.ToString(), [Text.UTF8Encoding]::new($false))
Write-Output "Generated $OutCpp (barrel=41 layers, endcap=39 layers, rgridb=$($rgridb.Count), zgride=$($zgride.Count))"
