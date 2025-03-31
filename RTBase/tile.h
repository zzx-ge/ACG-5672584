#pragma once

struct Tile {
    int xStart, xEnd;
    int yStart, yEnd;
    float variance;       // 当前瓦片的方差
    int currentSPP;       // 当前瓦片的采样数
    bool completed;       // 标记瓦片是否达到采样要求
};