##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/scripts/object.sh $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1917 $
## Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
## Description: shell script to update Makefile.objects files for templates
##

#
# Generate object information in the template subdrectories
#

PERL=${PERL:-perl}
DIFF="$PERL ../../scripts/cmp.pl"
OBJECTPL="$PERL ../../../scripts/object.pl"

SRCDIR=../$1
DIR=$2
rm -f $DIR/Makefile.objects.tmp

#
# Generate a Makefile.objects.tmp file in the specified directory
#

echo "Checking objects in directory "$DIR



(cd $DIR; $OBJECTPL libdefault \
    `cat $SRCDIR/default.filenames` \
    `cat $SRCDIR/Double.filenames` \
    `cat $SRCDIR/Integer.filenames`\
    `cat $SRCDIR/double.filenames` \
    `cat $SRCDIR/int.filenames` \
    )

(cd $DIR; $OBJECTPL libchar \
    `cat $SRCDIR/char.filenames` \
    )

(cd $DIR; $OBJECTPL libbool \
    `cat $SRCDIR/bool.filenames` \
    )

(cd $DIR; $OBJECTPL libfloat\
    `cat $SRCDIR/float.filenames` \
    `cat $SRCDIR/Float.filenames` \
    )
(cd $DIR; $OBJECTPL libdcomplex \
    `cat $SRCDIR/dcomplex.filenames` \
    `cat $SRCDIR/Complex.filenames` \
    )

#
# If Makefile.objects does not exist, then create it
#

if [ ! -r $DIR/Makefile.objects ] ; then
    echo "   creating "$DIR/Makefile.objects
    mv -f $DIR/Makefile.objects.tmp $DIR/Makefile.objects
    
#
# Otherwise, copy if the two files are not the same.  Remove the CVS
# portions of the header to ignore changes in date/revision/modified.
#

else
    if $DIFF $DIR/Makefile.objects.tmp $DIR/Makefile.objects; then
	rm -f $DIR/Makefile.objects.tmp
    else
	echo "   updating "$DIR/Makefile.objects
	mv -f $DIR/Makefile.objects.tmp $DIR/Makefile.objects
    fi
fi

exit 0
