#! /usr/bin/perl 

$usage = "Installer.pl -m hmmerDIR -p ProtexcluderDIR\n";

# to install Protexcluder

use Getopt::Std;

getopts("m:p:");

$hmmerDIR = defined $opt_m ? $opt_m : "";

$prexcDIR = defined $opt_p ? $opt_p : "";

@Raw_Files = glob "*.npl";
foreach(@Raw_Files) {
        $NPL = $_;
        open(RF, "$NPL")||die"$!\n";
        $PL = $NPL;
        $PL =~ s/\.npl/\.pl/;
        open(PL, ">$PL")||die"$!\n";
        while(<RF>) {
                chomp;
                $Line = $_;

                $Line =~ s/_hmmer_/$hmmerDIR/;
                $Line =~ s/_prexc_/$prexcDIR/;

                print(PL "$Line\n");
        }
        close(RF);
        close(PL);

        system "chmod 755 $PL\n";
}

print "Install finished!\n";
print "If you input the wrong path, you can do it again with corrected paths\n";
print "--------------------------- Have a nice day! -------------------------\n\n\n";

