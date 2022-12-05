void Engine::UpdateVoltages(unsigned int startX, unsigned int numX)
{
	unsigned int pos0, pos1, pos2;
  const unsigned int end1 = numLines[1];
  const unsigned int end2 = numLines[2];
	bool shift0, shift1, shift2;

	pos0 = startX;
	//voltage updates
	for (unsigned int posX=0; posX<numX; ++posX)
	{
		shift0=pos0;
		for (pos1=0; pos1<end1; ++pos1)
		{
			shift1=pos1;
      #pragma GCC unroll 4
			for (pos2=0; pos2<end2; ++pos2)
			{
				shift2=pos2;
				//do the updates here
				//for x
				volt[0][pos0][pos1][pos2] *= Op->vv[0][pos0][pos1][pos2];
				volt[0][pos0][pos1][pos2] += Op->vi[0][pos0][pos1][pos2] * ( curr[2][pos0][pos1][pos2] - curr[2][pos0][pos1-shift1][pos2] - curr[1][pos0][pos1][pos2] + curr[1][pos0][pos1][pos2-shift2]);

				//for y
				volt[1][pos0][pos1][pos2] *= Op->vv[1][pos0][pos1][pos2];
				volt[1][pos0][pos1][pos2] += Op->vi[1][pos0][pos1][pos2] * ( curr[0][pos0][pos1][pos2] - curr[0][pos0][pos1][pos2-shift2] - curr[2][pos0][pos1][pos2] + curr[2][pos0-shift0][pos1][pos2]);

				//for z
				volt[2][pos0][pos1][pos2] *= Op->vv[2][pos0][pos1][pos2];
				volt[2][pos0][pos1][pos2] += Op->vi[2][pos0][pos1][pos2] * ( curr[1][pos0][pos1][pos2] - curr[1][pos0-shift0][pos1][pos2] - curr[0][pos0][pos1][pos2] + curr[0][pos0][pos1-shift1][pos2]);
			}
		}
		++pos0;
	}
}
