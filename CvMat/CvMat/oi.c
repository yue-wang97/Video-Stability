#include <stdlib.h>  
#include<stdio.h>

int main()
{

FILE *pFile = fopen("test.txt", "r"); //��ȡ�ļ���ָ��
char *pBuf;  //�����ļ�ָ��
fseek(pFile, 0, SEEK_END); //��ָ���ƶ����ļ��Ľ�β ����ȡ�ļ�����
int len = ftell(pFile); //��ȡ�ļ�����
//pBuf =new char[len + 1]; //�������鳤��
pBuf=(char *)malloc(sizeof(char)*(len + 1));
rewind(pFile); //��ָ���ƶ����ļ���ͷ ��Ϊ����һ��ʼ��ָ���ƶ�����β��������ƶ����� �����
fread(pBuf, 1, len, pFile); //���ļ�
pBuf[len] = 0; //�Ѷ������ļ����һλ дΪ0 Ҫ��Ȼϵͳ��һֱѰ�ҵ�0��Ž���
fclose(pFile); // �ر��ļ�

FILE *fpt;
fpt = fopen("before50.txt", "w");//���ĵ���д�� 
fprintf(fpt, "%s", pBuf);
fclose(fpt);
}