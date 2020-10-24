#include <stdlib.h>  
#include<stdio.h>

int main()
{

FILE *pFile = fopen("test.txt", "r"); //获取文件的指针
char *pBuf;  //定义文件指针
fseek(pFile, 0, SEEK_END); //把指针移动到文件的结尾 ，获取文件长度
int len = ftell(pFile); //获取文件长度
//pBuf =new char[len + 1]; //定义数组长度
pBuf=(char *)malloc(sizeof(char)*(len + 1));
rewind(pFile); //把指针移动到文件开头 因为我们一开始把指针移动到结尾，如果不移动回来 会出错
fread(pBuf, 1, len, pFile); //读文件
pBuf[len] = 0; //把读到的文件最后一位 写为0 要不然系统会一直寻找到0后才结束
fclose(pFile); // 关闭文件

FILE *fpt;
fpt = fopen("before50.txt", "w");//打开文档，写入 
fprintf(fpt, "%s", pBuf);
fclose(fpt);
}