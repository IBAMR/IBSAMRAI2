#! /usr/bin/perl
##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/scripts/source_manipulation/findHeaderFiles.pl $
## Package:     SAMRAI scripts
## Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1917 $
## Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
## Description: perl script to update Xd sed files to templates on DIM
##

use File::Basename;
use File::Find;
use Cwd;

%directory_to_package = (
			 'apputils', 'appu',
			 'algorithm', 'algs',
			 'solvers', 'solv', 
			 'geometry', 'geom',
			 'multiblock', 'mblk',
			 'mesh', 'mesh',
			 'mathops', 'math',
			 'patchdata', 'pdat',
			 'transfer', 'xfer', 
			 'hierarchy', 'hier',
			 'toolbox', 'tbox', 
			 );

my $pwd = cwd;
#die basename($0) . " should not be run from your current directory"
#    if $pwd =~ m:\b(examples|source/test|source/scripts)(/|$):;


my $debug = 0;

#
# Read in sed.data to get the substitution patterns for X strings.
#

#
# Get prefix pattern to look for, or use default prefix pattern.
#

my $prepat;
$prepat = q|(.*\.h(.sed)?$)|;
print "prepat: $prepat\n" if ( $debug > 0 );


#
# Find the h files to search
#

@allfiles = ();
sub selectXfile {
    # This subroutine selects the X files in a find command.
    # print "-$File::Find::dir- -$File::Find::name-\n";
#    if ( $File::Find::dir =~ m!/(examples|source/test|source/scripts|CVS|[123]d)$!o ) {

    if ( $File::Find::dir =~ m!/(test|.svn|CVS|[123]d)$!o ) {
	# print "pruned\n";
	$File::Find::prune = true;
    }
    elsif ( -f && m/$prepat/o ) {
	push @allfiles, $File::Find::name;
	$allfiles[$#allfiles] =~ s|^\./||o;
    }
}
print "Scanning...\n" if ( $debug > 0 );
find( \&selectXfile, '.' );
print "Done.\n" if ( $debug > 0 );

my $debug=1;

for $f (@allfiles) {

    ($name,$path,$suffix) = fileparse($f,@suffixlist);

    ($packdir = $path) =~ s/(.*?)\/.*/\1/;

    
    $pack = $directory_to_package{$packdir};

    print "$pack $name\n" if $debug;

}
