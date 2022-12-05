void Engine::UpdateCurrents(unsigned int startX, unsigned int numX)
{
	unsigned int pos0, pos1, pos2;
  const unsigned int end1 = numLines[1]-1;
  const unsigned int end2 = numLines[2]-1;
	pos0 = startX;
	for (unsigned int posX=0; posX<numX; ++posX)
	{
		for (pos1=0; pos1<end1; ++pos1)
		{
      #pragma GCC unroll 4
			for (pos2=0; pos2<end2; ++pos2)
			{
				//do the updates here
				//for x
				curr[0][pos0][pos1][pos2] *= Op->ii[0][pos0][pos1][pos2];
				curr[0][pos0][pos1][pos2] += Op->iv[0][pos0][pos1][pos2] * ( volt[2][pos0][pos1][pos2] - volt[2][pos0][pos1+1][pos2] - volt[1][pos0][pos1][pos2] + volt[1][pos0][pos1][pos2+1]);

				//for y
				curr[1][pos0][pos1][pos2] *= Op->ii[1][pos0][pos1][pos2];
				curr[1][pos0][pos1][pos2] += Op->iv[1][pos0][pos1][pos2] * ( volt[0][pos0][pos1][pos2] - volt[0][pos0][pos1][pos2+1] - volt[2][pos0][pos1][pos2] + volt[2][pos0+1][pos1][pos2]);

				//for z
				curr[2][pos0][pos1][pos2] *= Op->ii[2][pos0][pos1][pos2];
				curr[2][pos0][pos1][pos2] += Op->iv[2][pos0][pos1][pos2] * ( volt[1][pos0][pos1][pos2] - volt[1][pos0+1][pos1][pos2] - volt[0][pos0][pos1][pos2] + volt[0][pos0][pos1+1][pos2]);
			}
		}
		++pos0;
	}
}
