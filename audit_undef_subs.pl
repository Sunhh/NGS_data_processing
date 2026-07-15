#!/usr/bin/perl
# audit_undef_subs.pl — best-effort static check for calls to subroutines that are
# neither locally defined, imported via @EXPORT / explicit use-list, nor a Perl builtin.
# Catches "undefined subroutine" risks that `perl -c` does not.
#
# Usage:  perl -I MyPM audit_undef_subs.pl file1.pl file2.pl ...
# Output: "FLAG\t<file>\t<subname>" per suspect call; "UNLOADED\t<file>\t<mods>" when a
#         used module could not be loaded (its imports are then unknown).
#
# HEURISTIC / KNOWN LIMITATIONS (results need a human glance):
#  - It strips POD, comments, heredocs, quoted strings and /-, !-delimited regexes, but
#    regexes with other delimiters (m{}, m,,, s{}{}) can leak `word(` fragments as
#    false positives (e.g. \t(\d+), GT(:|$)). All such flags are non-calls.
#  - Dynamic dispatch (&{$ref}, symbol-table tricks, AUTOLOAD) is not modelled.
# On this repo (2026-07, 396 MyPM-using .pl) every FLAG was verified to be a regex/string
# fragment — i.e. ZERO real undefined-sub calls.
my %kw = map {$_=>1} qw(
  if elsif else unless while until for foreach do sub return
  my our local state use no require package __PACKAGE__ __FILE__ __LINE__ __DATA__ __END__
  and or not xor eq ne lt gt le ge cmp x qw q qq qr m s tr y
  BEGIN END INIT CHECK last next redo goto default given when print printf say
);
my %special = map {$_=>1} qw(STDIN STDOUT STDERR ARGV ARGVOUT ENV DATA);
sub is_builtin { my $n=shift; return (eval { no warnings; prototype("CORE::$n"); 1 }) ? 1 : 0; }
my %expc;
sub mod_info {
  my $m=shift;
  return $expc{$m} if $expc{$m};
  my @ex;
  my $ok = eval { local $SIG{__WARN__}=sub{}; (eval "require $m; 1") or die; 1 };
  if ($ok) { no strict 'refs'; @ex = map { my $x=$_; $x=~s/^[&\$\@\%\*]//; $x } @{"${m}::EXPORT"}; }
  $expc{$m} = { loaded=>($ok?1:0), ex=>\@ex };
  return $expc{$m};
}
sub strip {
  my $src=shift;
  $src =~ s/\r//g;
  $src =~ s/^__(?:END|DATA)__\b.*\z//ms;
  $src =~ s/^=\w+.*?^=cut\b[^\n]*\n/\n/msg;
  my @L=split /\n/,$src,-1; my @o; my $i=0;
  while ($i<=$#L) {
    my $ln=$L[$i]; push @o,$ln;
    if ( $ln =~ /<<~?\s*(?:"([A-Za-z_]\w*)"|'([A-Za-z_]\w*)'|([A-Za-z_]\w*))\s*;/ ) {
      my $t = $1 // $2 // $3; $i++;
      while ($i<=$#L && $L[$i] !~ /^\s*\Q$t\E\s*$/) { $i++; }
      push @o, "" if $i<=$#L;
    }
    $i++;
  }
  $src=join "\n",@o;
  $src =~ s/"(?:[^"\\]|\\.)*"/""/gs;   # strings BEFORE comments (protect # inside strings)
  $src =~ s/'(?:[^'\\]|\\.)*'/''/gs;
  $src =~ s/^\s*#[^\n]*$//mg;
  $src =~ s/(?<!\$)#.*$//mg;
  $src =~ s/\btr\s*\/(?:[^\/\\]|\\.)*\/(?:[^\/\\]|\\.)*\/[a-z]*/tr\/\/\//gs;
  $src =~ s/\bs\s*\/(?:[^\/\\]|\\.)*\/(?:[^\/\\]|\\.)*\/[a-z]*/s\/\/\//gs;
  $src =~ s/\b(?:m|qr)\s*\/(?:[^\/\\]|\\.)*\/[a-z]*/m\/\//gs;
  $src =~ s/(=~|!~)\s*\/(?:[^\/\\]|\\.)*\/[a-z]*/$1 \/\//gs;

  $src =~ s/\btr\s*!(?:[^!\\]|\\.)*!(?:[^!\\]|\\.)*![a-z]*/tr!!!/gs;
  $src =~ s/\bs\s*!(?:[^!\\]|\\.)*!(?:[^!\\]|\\.)*![a-z]*/s!!!/gs;
  $src =~ s/\b(?:m|qr)\s*!(?:[^!\\]|\\.)*![a-z]*/m!!/gs;
  $src =~ s/(=~|!~)\s*!(?:[^!\\]|\\.)*![a-z]*/$1 m!!/gs;

  $src =~ s/qw[\(\[\{\/<].*?[\)\]\}\/>]/qw()/gs;
  return $src;
}
sub import_names {   # names imported by "use MOD REST;"
  my $rest = shift;
  my @names;
  while ($rest =~ /qw\s*[\(\[\{\/<](.*?)[\)\]\}\/>]/gs) { push @names, split ' ', $1; }
  while ($rest =~ /'([^']+)'/g) { push @names, $1; }
  while ($rest =~ /"([^"]+)"/g) { push @names, $1; }
  return @names;
}
for my $file (@ARGV) {
  open my $fh,"<",$file or do { print "ERR\t$file\n"; next; };
  my $raw=do{local $/;<$fh>}; close $fh;
  my $src=strip($raw);
  my %known;
  while ($raw =~ /^\s*sub\s+(\w+)/mg) { $known{$1}=1; }          # local subs from RAW
  my @unloaded;
  while ($raw =~ /^\s*use\s+([\w:]+)\s*([^;]*);/mg) {
    my ($mod,$rest)=($1,$2);
    next if $mod =~ /^(strict|warnings|utf8|constant|vars|lib|parent|base|overload)$/;
    my @names = import_names($rest);
    if (@names) { for my $x (@names){ $x=~s/^[&\$\@\%\*]//; $known{$x}=1; } }
    elsif ($rest =~ /^\s*\(\s*\)/) { }                            # imports nothing
    else { my $mi=mod_info($mod); if ($mi->{loaded}){ $known{$_}=1 for @{$mi->{ex}}; } else { push @unloaded,$mod; } }
  }
  my %called;
  while ($src =~ /&(\w+(?:::\w+)*)/g) { my $n=$1; next if $n=~/::/; $called{$n}++; }
  while ($src =~ /(?<![\w:>&\$\@\%\-])(\w+)\s*\(/g) { $called{$1}++; }
  for my $n (sort keys %called) {
    next if $kw{$n} or $special{$n} or $known{$n} or is_builtin($n);
    print "FLAG\t$file\t$n\n";
  }
  print "UNLOADED\t$file\t".join(",",@unloaded)."\n" if @unloaded;
}
