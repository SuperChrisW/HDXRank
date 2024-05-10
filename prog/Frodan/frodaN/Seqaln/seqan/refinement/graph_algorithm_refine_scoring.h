 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: graph_algorithm_refine_scoring.h,v 1.1 2011/04/23 04:46:55 dwfarrel Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_REFINE_SCORING_H
#define SEQAN_HEADER_GRAPH_REFINE_SCORING_H


namespace SEQAN_NAMESPACE_MAIN
{

	



//fake score function 
template<typename TScoreValue,typename TStringSet,typename TAlign,typename TValue, typename TSize>
TScoreValue
getScore(TScoreValue &,
		 TStringSet &,
		 TAlign &,
		 TValue,
		 TValue,
		 TSize,
		 TSize)
{
SEQAN_CHECKPOINT
	return 1;
}				




}
#endif //#ifndef SEQAN_HEADER_...