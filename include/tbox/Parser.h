//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/inputdb/Parser.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Parser that reads the input database grammar
//

#ifndef included_tbox_Parser
#define included_tbox_Parser

#include "SAMRAI_config.h"
#ifndef included_stdio
#define included_stdio
#include <stdio.h>
#endif
#include "tbox/Database.h"
#include "tbox/List.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
#define included_String
#endif

namespace SAMRAI {
   namespace tbox {

/**
 * Class Parser parses the user input file and places the resulting
 * (key,value) pairs into the input database object.  The parser object
 * controls the overall parsing of the input file, which includes error
 * handing and tracking file, line number, and cursor position.  If
 * running on multiple processors, only node zero reads in data from the
 * specified input file and broadcasts that data to the other processors.
 * The input file argument for the other processors is ignored and may be
 * NULL.
 *
 * The parser class also defines a ``default'' parser that may be accessed
 * via a static member function.  The default parser may only be accessed
 * by the yacc/lex routines called during the input file parsing.  This
 * singleton-like approach provides a clean way to communicate parser
 * information to the yacc/lex files without global variables.
 *
 * This parser (and the associated yacc and lex files) assume the GNU flex
 * and bison utilities.  These utilities are only required if the grammar or
 * scanner is changed, since the SAMRAI distribution already contains the
 * output files from flex and bison.
 */

class Parser
{
public:
   /**
    * The parser constructor simply creates an uninitialized parser object.
    * Member function parse() must be called before any other member function
    * to initialize the object and parse the input data.  Function parse()
    * may be called multiple times to parse multiple input files, but all
    * state values (such as the number of errors or warnings) are reset at
    * the beginning of each new parse pass.
    */
   Parser();

   /**
    * Destroy the parser object and deallocate parser data.
    */
   ~Parser();

   /**
    * Parse the input file from the specified file stream.  The number of
    * syntax errors is returned.  A successful parse will return zero errors.
    * The parse() function takes the initial filename (for informational
    * purposes) and the filestream from which to read the parse data.
    * All (key,value) pairs are placed in the specified database.  If
    * running in parallel, the fstream must be valid on node zero, but
    * is ignored on other nodes and may be set to NULL.  Multiple input
    * files may be parsed by calling parse() for each file, but all variables
    * are reset at the beginning of each parse.
    */
   int parse(
      const std::string& filename,
      FILE* fstream,
      Pointer<Database> database);

   /**
    * Return the total number of errors resulting from the parse.
    */
   int getNumberErrors() const;

   /**
    * Return the total number of warnings resulting from the parse.
    */
   int getNumberWarnings() const;

   /**
    * Return the parser object.  This mechanism is useful for communicating
    * with the yacc/lex routines during the input file parse.  The default
    * parser will be NULL outside of the parse call.
    */
   static Parser *getParser();

   /**
    * Return the current database scope.  The current scope is modified
    * through the enterScope() and leaveScope() member functions.
    */
   Pointer<Database>& getScope();

   /**
    * Create a new database scope with the specified name.  This new scope
    * will be the default scope until leaveScope() is called.
    */
   void enterScope(const std::string& name);

   /**
    * Leave the current database scope and return to the previous scope.
    * It is an error to leave the outermost scope.
    */
   void leaveScope();

   /**
    * Lookup the scope that contains the specified key.  If the scope does
    * not exist, then return a NULL pointer to the database.
    */
   Pointer<Database> getDatabaseWithKey(const std::string& name);

   /**
    * Save the current context and switch to the specified input file.
    * This routine returns true if the file exists and the switch was
    * successful and false otherwise.
    */
   bool pushIncludeFile(const std::string& filename);

   /**
    * Pop the include file context off of the stack and return to the
    * previous include file.
    */
   void popIncludeFile();

   /**
    * Report a parsing error with the specified error message.  This routine
    * will only be called from the parser or the scanner.  Errors are printed
    * to pout, since it is assumed that all nodes are parsing the same input
    * file.
    */
   void error(const std::string& message);

   /**
    * Report a parsing warning with the specified warning message.  This
    * routine will only be called from the parser or the scanner.  Errors
    * are printed to pout, since it is assumed that all nodes are parsing
    * the same input file.
    */
   void warning(const std::string& message);

   /**
    * Set the input line which is currently being parsed.
    */
   void setLine(const std::string& line);

   /**
    * Advance the line number by the specified number of lines.  If no
    * argument is given, then the line number is advanced by one line.
    */
   void advanceLine(const int nline = 1);

   /**
    * Advance the position of the cursor on the line using the values
    * in the specified character string.  Tab characters in the string
    * are assumed to advance the cursor to eight character tab stops.
    * The cursor position is automatically reset to one whenever the
    * line number is changed.
    */
   void advanceCursor(const std::string& token);

   /**
    * Define the input reading routine used by flex.  Under MPI, node zero
    * reads the input and broadcasts the character data to all processors.
    */
   int yyinput(char *buffer, const int max_size);

private:
   Parser(const Parser&);	// not implemented
   void operator=(const Parser&);	// not implemented

   struct ParseData {
      std::string d_filename;	// filename for description
      FILE* d_fstream;		// input stream to parse
      std::string d_linebuffer;      // line being parsed
      int d_linenumber;		// line number in input stream
      int d_cursor;		// cursor position in line
      int d_nextcursor;		// next cursor position in line
   };

   int d_errors;		// total number of parse errors
   int d_warnings;		// total number of warnings

#ifdef LACKS_NAMESPACE_IN_DECLARE
   List< ParseData > d_parse_stack;
#else
   List< Parser::ParseData > d_parse_stack;
#endif

   List< Pointer< Database > > d_scope_stack;

   static Parser* s_default_parser;
   
   static bool s_static_tables_initialized;

   std::string d_pathname;           // path to filename for including
};


}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/Parser.I"
#endif
#endif
