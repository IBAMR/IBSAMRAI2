%{
//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-4/source/toolbox/inputdb/Scanner.l $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Lex scanner description for the SAMRAI input database
//

// Ignore Intel warnings about unreachable statements
#ifdef __INTEL_COMPILER
#pragma warning (disable:177)
#endif

#include "SAMRAI_config.h"
#include <stdlib.h>
#include "tbox/Complex.h"
#include "tbox/Parser.h"
#include <string>
using namespace std;

#define YY_NEVER_INTERACTIVE 1
#define YY_NO_UNPUT 1

using namespace SAMRAI;
using namespace tbox;

struct KeyData;
#include "Grammar.h"

#define YY_INPUT(b,r,s) r = Parser::getParser()->yyinput(b, s)

#define TOKEN(T) {						\
   Parser::getParser()->advanceCursor(yytext);	\
   return(T);							\
}

#define MAXIMUM_INCLUDE_DEPTH 16
static YY_BUFFER_STATE s_include_stack[MAXIMUM_INCLUDE_DEPTH];
static int s_include_stack_top = 0;

%}

%option noyywrap
%x COMMENT

S_CHAR    (\'[^\']\')
S_INTEGER ([-]?[0-9]+)
S_KEYWORD ([a-zA-Z][a-zA-Z0-9_]*)
S_DOUBLE  (([-]?[0-9]+"."[0-9]*([Ee][+-]?[0-9]+)?)|([-]?"."[0-9]*([Ee][+-]?[0-9]+)?)|([-]?[0-9]+([Ee][+-]?[0-9]+)))
S_STRING  (\"([^\n"])*\")

%%

"&&"    TOKEN(T_AND) 
"="     TOKEN(T_ASSIGN)
"}"     TOKEN(T_CLOSE_CURLY)
")"     TOKEN(T_CLOSE_PAREN)
"]"     TOKEN(T_CLOSE_SQBKT)
","     TOKEN(T_COMMA)
"/"     TOKEN(T_DIV) 
":"     TOKEN(T_ELSE)
"=="    TOKEN(T_EQUALS)
"^"     TOKEN(T_EXP) 
"FALSE" TOKEN(T_FALSE)
">="    TOKEN(T_GREATER_EQUALS)
">"     TOKEN(T_GREATER)
"<="    TOKEN(T_LESS_EQUALS)
"<"     TOKEN(T_LESS)
"-"     TOKEN(T_MINUS) 
"*"     TOKEN(T_MULT) 
"!"     TOKEN(T_NOT) 
"!="    TOKEN(T_NOT_EQUALS)
"||"    TOKEN(T_OR) 
"{"     TOKEN(T_OPEN_CURLY) 
"("     TOKEN(T_OPEN_PAREN) 
"["     TOKEN(T_OPEN_SQBKT)
"+"     TOKEN(T_PLUS) 
"?"     TOKEN(T_QUESTION)
";"     TOKEN(T_SEMI)
"TRUE"  TOKEN(T_TRUE)

{S_CHAR} {
   Parser::getParser()->advanceCursor(yytext);
   yylval.u_char = yytext[1];
   return(T_CHAR);
}

{S_INTEGER} {
   Parser::getParser()->advanceCursor(yytext);
   yylval.u_integer = atol(yytext);
   return(T_INTEGER);
}

{S_KEYWORD} {
   Parser::getParser()->advanceCursor(yytext);
   yylval.u_keyword = new string(yytext);
   return(T_KEYWORD);
}

{S_DOUBLE} {
   Parser::getParser()->advanceCursor(yytext);
   yylval.u_double = atof(yytext);
   return(T_DOUBLE);
}

{S_STRING} {
   Parser::getParser()->advanceCursor(yytext);
   string s(yytext);
   yylval.u_string = new string(s.substr(1,s.length()-2));
   return(T_STRING);
}

[ \t]+ {
   Parser::getParser()->advanceCursor(yytext);
}

\n.* {
   Parser::getParser()->setLine(yytext+1);
   Parser::getParser()->advanceLine();
   yyless(1);
}

"#include"[ \t]+\"[^\"\n]+\"[ \t]*\n {
   string s(yytext);
   string::size_type start = s.find("\"")+1;
   string::size_type end   = s.rfind("\"")-1;

   string filename(s, start, end-start+1);

   if (s_include_stack_top >= MAXIMUM_INCLUDE_DEPTH) {
      Parser::getParser()->error("Too many nested #include statements");
      Parser::getParser()->advanceLine();
   } else {
      if (Parser::getParser()->pushIncludeFile(filename)) {
         s_include_stack[s_include_stack_top++] = YY_CURRENT_BUFFER;
         yy_switch_to_buffer(yy_create_buffer(NULL,YY_BUF_SIZE));
      }
   }
}

"//".*\n {
   Parser::getParser()->advanceLine();
}

"/*" {
   Parser::getParser()->advanceCursor(yytext);
   BEGIN COMMENT;
}

<COMMENT>"/*" {
   Parser::getParser()->advanceCursor(yytext);
   Parser::getParser()->warning("``/*'' found within comment");
}

<COMMENT>\n {
   Parser::getParser()->advanceLine();
}

<COMMENT>"*/" {
   Parser::getParser()->advanceCursor(yytext);
   BEGIN INITIAL;
}

<COMMENT><<EOF>> {
   Parser::getParser()->warning("EOF detected within comment");
   yy_delete_buffer(yy_current_buffer);
   yyterminate();
}

<COMMENT>. {
   Parser::getParser()->advanceCursor(yytext);
}

. {
   Parser::getParser()->advanceCursor(yytext); 
   Parser::getParser()->warning("Illegal character token in input");
   return(yytext[0]);
}

<<EOF>> {
   yy_delete_buffer(yy_current_buffer);
   if (s_include_stack_top == 0) {
      yyterminate();
   } else {
      yy_switch_to_buffer(s_include_stack[--s_include_stack_top]);
      Parser::getParser()->popIncludeFile();
      Parser::getParser()->advanceLine();
   }
}

%%
