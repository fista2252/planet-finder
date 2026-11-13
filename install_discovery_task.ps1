param(
    [Parameter(Mandatory = $true)][string]$Python,
    [Parameter(Mandatory = $true)][string]$Config
)

$ErrorActionPreference = 'Stop'
$pythonPath = (Resolve-Path -LiteralPath $Python).Path
$configPath = (Resolve-Path -LiteralPath $Config).Path
$worker = Join-Path $PSScriptRoot 'discovery_worker.py'
$name = 'Planet Finder Discovery'
if ((Get-TimeZone).Id -ne 'Pakistan Standard Time') {
    throw 'This installer expects the Windows time zone Pakistan Standard Time.'
}
if (Get-ScheduledTask -TaskName $name -ErrorAction SilentlyContinue) {
    throw 'The discovery task already exists. Inspect its configuration before replacing it.'
}
$action = New-ScheduledTaskAction -Execute $pythonPath -Argument ('"{0}" tick --config "{1}"' -f $worker, $configPath) -WorkingDirectory $PSScriptRoot
$trigger = New-ScheduledTaskTrigger -Weekly -DaysOfWeek Tuesday,Wednesday,Thursday,Friday,Saturday,Sunday -At '16:00'
$settings = New-ScheduledTaskSettingsSet -MultipleInstances IgnoreNew -ExecutionTimeLimit (New-TimeSpan -Minutes 10)
$principal = New-ScheduledTaskPrincipal -UserId ([Security.Principal.WindowsIdentity]::GetCurrent().Name) -LogonType Interactive
Register-ScheduledTask -TaskName $name -Action $action -Trigger $trigger -Settings $settings -Principal $principal -Description 'Deterministic Planet Finder archive work. No model invocation.'
