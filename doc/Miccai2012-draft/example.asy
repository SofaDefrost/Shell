pair quadSize = (100, 100);
pen defaultPen = black+beveljoin+1;

void drawQuad(pair pos, pair size, pen style = defaultPen)
{
    draw((pos.x, pos.y)--(pos.x + size.x, pos.y)--(pos.x + size.x, pos.y + size.y)--(pos.x, pos.y + size.y)--(pos.x, pos.y)--(pos.x + size.x, pos.y + size.y), style);
}

void drawQuadArray(pair pos, pair size, int width, int height, pen style = defaultPen)
{
    for (int x = 0; x < width; ++x)
    {
        for (int y = 0; y < height; ++y)
        {
            drawQuad((pos.x + x * size.x, pos.y + y * size.y), size, style);
        }
    }
}

drawQuadArray((0, 0), quadSize, 3, 3);
drawQuadArray((quadSize.x * 5, 0), quadSize, 3, 3);

for (int y = 0; y < 4; ++y)
{
    dot((quadSize.x * 3, y * quadSize.y), blue+7);
    draw((quadSize.x * 3, y * quadSize.y)--(quadSize.x * 5, y * quadSize.y), blue+dashed+beveljoin+1, Arrow);
}

drawQuadArray((0, quadSize.y * 4), quadSize, 2, 3);
drawQuadArray((quadSize.x * 5, quadSize.y * 4), quadSize, 3, 3);
drawQuadArray((quadSize.x * 2, quadSize.y * 4), (quadSize.x * 3, quadSize.y), 1, 3, blue+beveljoin+5);