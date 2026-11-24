# SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
# SPDX-License-Identifier: GPL-3.0-or-later
<#
.SYNOPSIS
  Zero-dependency PowerShell companion to tools/analyze_simana.py.
  Generates inventory.json, commons.json, banks.json, and SUMMARY.md from simana_v94d/src/car.
#>
param(
  [Parameter(Mandatory = $true)][string]$CarDir,
  [Parameter(Mandatory = $true)][string]$OutDir
)

$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

$rePatch   = [regex]::new('^\+PATCH\s*,\s*([^.,\s]+)', 'IgnoreCase,Compiled')
$reDeck    = [regex]::new('^\+DECK\s*,\s*([^.,\s]+)', 'IgnoreCase,Compiled')
$reKeep    = [regex]::new('^\+KEEP\s*,\s*([^.,\s]+)', 'IgnoreCase,Compiled')
$reCde     = [regex]::new('^\+(?:CDE|SEQ)\s*,\s*([^.,\s]+)', 'IgnoreCase,Compiled')
$reSelfIf  = [regex]::new('^\+SELF\s*,.*IF\s*=', 'IgnoreCase,Compiled')
$reRoutine = [regex]::new('^\s{6}\s*(?:SUBROUTINE|(?:[A-Z0-9*]+\s+)?FUNCTION|PROGRAM|BLOCK\s*DATA)\s+([A-Z0-9_]+)', 'IgnoreCase,Compiled')
$reCommon  = [regex]::new('\bCOMMON\s*/\s*([A-Z0-9_]+)\s*/', 'IgnoreCase,Compiled')
$reMz      = [regex]::new('\bCALL\s+(MZBOOK|MZLIFT|MZDROP|MZGARB|MZWIPE|MZPUSH)\b', 'IgnoreCase,Compiled')
$reZebra   = [regex]::new('\b(LQ|IQ|Q)\s*\(', 'IgnoreCase,Compiled')
$reGoto    = [regex]::new('\bGO\s*TO\b', 'IgnoreCase,Compiled')
$reEquiv   = [regex]::new('\bEQUIVALENCE\b', 'IgnoreCase,Compiled')
$reEntry   = [regex]::new('^\s{6}\s*ENTRY\s+([A-Z0-9_]+)', 'IgnoreCase,Compiled')

$carFiles = Get-ChildItem -Path $CarDir -Filter '*.car' | Sort-Object Name
$inventory = New-Object System.Collections.Generic.List[object]
$commonsMap = @{}
$keepsMap = @{}
$mzCounts = @{}

foreach ($f in $carFiles) {
  $lines = [IO.File]::ReadAllLines($f.FullName, [Text.Encoding]::GetEncoding(28591))
  $currentPatch = ''
  $currentKeep = ''
  $codeLines = 0
  $patchyDir = 0
  $patches = 0
  $decks = 0
  $keeps = 0
  $cdeInc = 0
  $selfIf = 0
  $entries = 0
  $zebraRefs = 0
  $mzCalls = 0
  $gotos = 0
  $equivs = 0
  $routines = New-Object System.Collections.Generic.List[string]
  $fileCommons = New-Object System.Collections.Generic.HashSet[string]

  for ($i = 0; $i -lt $lines.Length; $i++) {
    $line = $lines[$i]
    if ([string]::IsNullOrWhiteSpace($line)) { continue }
    $ch = $line[0]
    if ($ch -eq '+') {
      $patchyDir++
      $m = $rePatch.Match($line)
      if ($m.Success) { $currentPatch = $m.Groups[1].Value.ToUpper(); $currentKeep = ''; $patches++; continue }
      $m = $reDeck.Match($line)
      if ($m.Success) { $currentKeep = ''; $decks++; continue }
      $m = $reKeep.Match($line)
      if ($m.Success) {
        $currentKeep = $m.Groups[1].Value.ToUpper()
        $keeps++
        if (-not $keepsMap.ContainsKey($currentKeep)) {
          $keepsMap[$currentKeep] = @{
            defined_in = $f.Name; patch = $currentPatch; line = ($i + 1)
            commons = (New-Object System.Collections.Generic.HashSet[string])
            used_in = (New-Object System.Collections.Generic.HashSet[string])
          }
        } elseif ($null -eq $keepsMap[$currentKeep].defined_in) {
          $keepsMap[$currentKeep].defined_in = $f.Name
          $keepsMap[$currentKeep].patch = $currentPatch
          $keepsMap[$currentKeep].line = ($i + 1)
        }
        continue
      }
      $m = $reCde.Match($line)
      if ($m.Success) {
        $cdeInc++
        $seq = $m.Groups[1].Value.ToUpper()
        if (-not $keepsMap.ContainsKey($seq)) {
          $keepsMap[$seq] = @{
            defined_in = $null; patch = ''; line = 0
            commons = (New-Object System.Collections.Generic.HashSet[string])
            used_in = (New-Object System.Collections.Generic.HashSet[string])
          }
        }
        [void]$keepsMap[$seq].used_in.Add($f.Name)
      }
      if ($reSelfIf.IsMatch($line)) { $selfIf++ }
      continue
    }
    if ($ch -eq 'C' -or $ch -eq 'c' -or $ch -eq '*' -or $ch -eq '!') { continue }

    $codeLines++
    $m = $reRoutine.Match($line)
    if ($m.Success) { $routines.Add($m.Groups[1].Value.ToUpper()) }
    elseif ($reEntry.IsMatch($line)) { $entries++ }

    foreach ($cm in $reCommon.Matches($line)) {
      $cn = $cm.Groups[1].Value.ToUpper()
      [void]$fileCommons.Add($cn)
      if (-not $commonsMap.ContainsKey($cn)) {
        $commonsMap[$cn] = @{
          files = (New-Object System.Collections.Generic.HashSet[string])
          keeps = (New-Object System.Collections.Generic.HashSet[string])
        }
      }
      [void]$commonsMap[$cn].files.Add($f.Name)
      if ($currentKeep -ne '') {
        [void]$commonsMap[$cn].keeps.Add($currentKeep)
        [void]$keepsMap[$currentKeep].commons.Add($cn)
      }
    }

    $mzm = $reMz.Match($line)
    if ($mzm.Success) {
      $mzCalls++
      $mzName = $mzm.Groups[1].Value.ToUpper()
      if (-not $mzCounts.ContainsKey($mzName)) { $mzCounts[$mzName] = 0 }
      $mzCounts[$mzName]++
    }
    $zebraRefs += $reZebra.Matches($line).Count
    if ($reGoto.IsMatch($line)) { $gotos++ }
    if ($reEquiv.IsMatch($line)) { $equivs++ }
  }

  $sortedCommons = @($fileCommons | Sort-Object)
  $inventory.Add([pscustomobject][ordered]@{
    file = $f.Name
    stem = $f.BaseName.ToLower()
    total_lines = $lines.Length
    code_lines = $codeLines
    patchy_directives = $patchyDir
    patches = $patches
    decks = $decks
    keeps = $keeps
    cde_includes = $cdeInc
    self_if_guards = $selfIf
    routine_count = $routines.Count
    entries = $entries
    common_count = $sortedCommons.Count
    commons = $sortedCommons
    zebra_array_refs = $zebraRefs
    mz_calls = $mzCalls
    gotos = $gotos
    equivalences = $equivs
  })
}

$commonsOut = [ordered]@{}
foreach ($k in ($commonsMap.Keys | Sort-Object)) {
  $commonsOut[$k] = [ordered]@{
    file_count = $commonsMap[$k].files.Count
    files = @($commonsMap[$k].files | Sort-Object)
    keeps = @($commonsMap[$k].keeps | Sort-Object)
  }
}

$keepsOut = [ordered]@{}
foreach ($k in ($keepsMap.Keys | Sort-Object)) {
  $keepsOut[$k] = [ordered]@{
    defined_in = $keepsMap[$k].defined_in
    patch = $keepsMap[$k].patch
    line = $keepsMap[$k].line
    commons = @($keepsMap[$k].commons | Sort-Object)
    used_in_file_count = $keepsMap[$k].used_in.Count
  }
}

$inventory | ConvertTo-Json -Depth 4 | Set-Content -Encoding UTF8 (Join-Path $OutDir 'inventory.json')
$commonsOut | ConvertTo-Json -Depth 4 | Set-Content -Encoding UTF8 (Join-Path $OutDir 'commons.json')
([ordered]@{ keeps = $keepsOut; mz_call_counts = $mzCounts }) | ConvertTo-Json -Depth 4 | Set-Content -Encoding UTF8 (Join-Path $OutDir 'banks.json')

$totLines = ($inventory | Measure-Object -Property total_lines -Sum).Sum
$totCode  = ($inventory | Measure-Object -Property code_lines -Sum).Sum
$totRout  = ($inventory | Measure-Object -Property routine_count -Sum).Sum

$detAna = $inventory | Where-Object { $_.stem -like '*ana' -and $_.stem -notin @('delana43','dstana10','hctanaxx') } | Sort-Object code_lines, common_count
$topCommons = $commonsOut.GetEnumerator() | Sort-Object { $_.Value.file_count } -Descending | Select-Object -First 20

$md = New-Object System.Collections.Generic.List[string]
$md.Add('# SIMANA v94D Static Inventory Summary')
$md.Add('')
$md.Add("- **``.car`` files:** $($inventory.Count)")
$md.Add("- **Total lines:** $('{0:N0}' -f $totLines) ($('{0:N0}' -f $totCode) non-comment code lines)")
$md.Add("- **Routines:** $('{0:N0}' -f $totRout)")
$md.Add("- **Distinct COMMON blocks:** $('{0:N0}' -f $commonsOut.Count)")
$md.Add("- **Distinct PATCHY ``+KEEP`` / ``+CDE`` sequences:** $('{0:N0}' -f $keepsOut.Count)")
$md.Add('')
$md.Add('## Phase 2 Pilot Candidates (Subdetector `*ana.car`, smallest first)')
$md.Add('')
$md.Add('| Package | Code lines | Routines | COMMONs | ZEBRA refs | GOTOs |')
$md.Add('|---|---:|---:|---:|---:|---:|')
foreach ($r in $detAna) {
  $md.Add("| ``$($r.file)`` | $('{0:N0}' -f $r.code_lines) | $($r.routine_count) | $($r.common_count) | $('{0:N0}' -f $r.zebra_array_refs) | $('{0:N0}' -f $r.gotos) |")
}
$md.Add('')
$md.Add('## Most Widely Shared COMMON Blocks (Top 20 by file count)')
$md.Add('')
$md.Add('| COMMON | Files | Defined in `+KEEP` |')
$md.Add('|---|---:|---|')
foreach ($kv in $topCommons) {
  $ks = if ($kv.Value.keeps.Count -gt 0) { ($kv.Value.keeps | Select-Object -First 4) -join ', ' } else { '*(inline)*' }
  $md.Add("| ``/$($kv.Key)/`` | $($kv.Value.file_count) | $ks |")
}

$md | Set-Content -Encoding UTF8 (Join-Path $OutDir 'SUMMARY.md')
Write-Output "Inventory written to $OutDir (files=$($inventory.Count), routines=$totRout, commons=$($commonsOut.Count), keeps=$($keepsOut.Count))"
