### Basic functions.
function exe_cmd {
  echo "[$(date)][CMD] $1"
  eval "$1"
  if [[ $? -eq "0" ]]
  then
    echo "[$(date)][CMD_done] $1"
  else
    echo "[$(date)][CMD_err] $1"
  fi
}

function tsmsg {
  echo "[$(date)][Msg] $1"
}


