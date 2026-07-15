#!/bin/bash
# Regenerate per-directory README.md script indexes + top-level SCRIPTS_INDEX.md.
# Each script's synopsis is pulled from its own `perl $0` / `Usage: $0` line.
#
# Auto-generated files carry the MARKER line below. Only files that are MISSING or
# carry the MARKER are (re)written; any README.md WITHOUT the marker is treated as
# hand-written and left untouched (and inlined into the master index for coverage).
# Idempotent and independent of git state. Usage:  bash gen_script_index.sh
set -u
cd "$(dirname "$0")"
DATE=$(date +%F)
MARKER='<!-- AUTOGEN gen_script_index.sh; do not edit — re-run the generator instead -->'

read -r -d '' EXTRACT <<'PL'
  use strict; use warnings;
  my $file=shift;
  open my $fh,"<",$file or die; my @L=<$fh>; close $fh; chomp @L;
  my $joined=join("\n",@L);
  my $usage='';
  for my $re ( qr/perl\s+\$0([^\n"']*)/, qr/[Uu]sage[:\s]+\$0([^\n"']*)/,
               qr/\bprint[^\n]*?\$0([^\n"']*)/, qr/\bdie[^\n]*?\$0([^\n"']*)/ ){
    if($joined =~ $re){
      my $u=$1; $u=~s/\\t/ /g; $u=~s/\\n.*$//; $u=~s/\s+/ /g; $u=~s/^\s+|\s+$//g;
      if(length($u)>=1){ $usage=$u; last; }
    }
  }
  my $desc='';
  if($usage ne ''){
    $desc="\$0 $usage";
  } else {
    for my $i (1..$#L){
      last if $L[$i] =~ /^\s*(use|my|our|package|require|GetOptions|BEGIN)\b/;
      if($L[$i] =~ /^#+\s*(.+?)\s*$/){
        my $c=$1;
        next if $c=~/^[-=*#\s]*$/;
        next if $c=~m{^#?\!} or $c=~m{/usr/bin|/bin/env};
        next if $c=~/^\d{4}[-\/.]\d/;
        next if $c=~m{https?://};
        next if $c=~/^(version|ver|copy|add(s|ed)?\s|Author|email|Date|Ref(erence)?|See|Note|TODO|FIXME|edit)\b/i;
        next if length($c)<8;
        $desc=$c; last;
      }
    }
  }
  $desc =~ s/\|/\\|/g;
  $desc = '_(no synopsis found)_' if $desc eq '';
  print "$desc\n";
PL

emit_table () {
  local d="$1" f b syn
  echo "| script | synopsis |"; echo "|---|---|"
  for f in $(ls "$d"/*.pl 2>/dev/null | sort); do
    b=$(basename "$f"); syn=$(perl -e "$EXTRACT" "$f"); syn=${syn//\$0/$b}
    printf '| `%s` | %s |\n' "$b" "$syn"
  done
}
anchor () { echo "$1" | sed 's|[/.]|-|g'; }
# hand-written = README.md exists AND does not carry our MARKER
is_handwritten () { [ -f "$1/README.md" ] && ! grep -qF "$MARKER" "$1/README.md"; }

mapfile -t DIRS < <(find . -type f -name '*.pl' | sed 's|/[^/]*$||' | sed 's|^\./||' | sort -u)
TOTP=$(find . -type f -name '*.pl' | wc -l); TOTD=${#DIRS[@]}

# 1) per-subdir READMEs (skip hand-written; always skip top-level ".")
for d in "${DIRS[@]}"; do
  [ "$d" = "." ] && continue
  is_handwritten "$d" && continue
  cnt=$(ls "$d"/*.pl 2>/dev/null | wc -l)
  { echo "$MARKER"
    echo "# \`$d\` — script index"; echo
    echo "_${cnt} Perl scripts. Auto-generated ${DATE}; synopsis is from each script's own \`perl \$0\` / \`Usage: \$0\` line — run a script with \`-h\` or no args for full options._"; echo
    emit_table "$d"; } > "$d/README.md"
done

# 2) master index — inline the dirs (with .pl) that are hand-written or top-level
INLINE="."; for d in "${DIRS[@]}"; do { [ "$d" != "." ] && is_handwritten "$d"; } && INLINE="$INLINE $d"; done
{
  echo "$MARKER"
  echo "# NGS_data_processing — master script index"; echo
  echo "_Auto-generated ${DATE}. **${TOTP} Perl scripts** across **${TOTD} directories**. Synopsis is from each script's own \`perl \$0\` / \`Usage: \$0\` line — run a script with \`-h\` or no args for full options._"; echo
  echo "## Directory overview"; echo
  echo "| directory | #scripts | index |"; echo "|---|---:|---|"
  for d in "${DIRS[@]}"; do
    cnt=$(ls "$d"/*.pl 2>/dev/null | wc -l)
    if [ "$d" = "." ] || is_handwritten "$d"; then
      lbl="README"; [ "$d" = "." ] && lbl="top README"
      printf '| `%s` | %d | [%s](%s/README.md) · [scripts](#%s) |\n' "$d" "$cnt" "$lbl" "$d" "$(anchor "$d")"
    else
      printf '| `%s` | %d | [README](%s/README.md) |\n' "$d" "$cnt" "$d"
    fi
  done
  echo; echo "> Directories with a hand-written README are not auto-overwritten; their script listings are inlined below for completeness."; echo
  echo "## Script listing — hand-written-README & top-level directories"
  for d in $INLINE; do
    cnt=$(ls "$d"/*.pl 2>/dev/null | wc -l); dn="$d"; [ "$d" = "." ] && dn="(top level)"
    echo; echo "<a id=\"$(anchor "$d")\"></a>"; echo "### \`$dn\` ($cnt)"; echo; emit_table "$d"
  done
} > SCRIPTS_INDEX.md
echo "done: ${TOTP} scripts / ${TOTD} dirs; inlined dirs: $INLINE"
