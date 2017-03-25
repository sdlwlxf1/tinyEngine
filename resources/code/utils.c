//
//  utils.h
//  tinyEngine
//
//  Created by lixuefeng on 2017/3/9.
//  Copyright © 2017年 sdlwlxf. All rights reserved.
//
#include "utils.h"
#include "png.h"
#include "tinyobj_loader_c.h"

#define FLT_MIN 1.175494351e-38F
#define FLT_MAX 3.402823466e+38F
#define O_RDONLY         00
#define O_WRONLY         01
#define O_RDWR           02

void CalcNormal(float N[3], float v0[3], float v1[3], float v2[3]) {
    float v10[3];
    float v20[3];
    float len2;
    
    v10[0] = v1[0] - v0[0];
    v10[1] = v1[1] - v0[1];
    v10[2] = v1[2] - v0[2];
    
    v20[0] = v2[0] - v0[0];
    v20[1] = v2[1] - v0[1];
    v20[2] = v2[2] - v0[2];
    
    N[0] = v20[1] * v10[2] - v20[2] * v10[1];
    N[1] = v20[2] * v10[0] - v20[0] * v10[2];
    N[2] = v20[0] * v10[1] - v20[1] * v10[0];
    
    len2 = N[0] * N[0] + N[1] * N[1] + N[2] * N[2];
    if (len2 > 0.0f) {
        float len = (float)sqrt((double)len2);
        
        N[0] /= len;
        N[1] /= len;
    }
}

int load_obj(float bmin[3], float bmax[3], const char* name, const char* type) {
    tinyobj_attrib_t attrib;
    tinyobj_shape_t* shapes = NULL;
    size_t num_shapes;
    tinyobj_material_t* materials = NULL;
    size_t num_materials;
    
    const char* ext = strrchr(type, '.');
    
    if (strcmp(ext, ".gz") == 0) {
        assert(0); /* todo */
    } else {
        //data = mmap_file(&data_len, filename);
        FILE * pFile;
        long lSize;
        size_t data_len;
        char * buffer;
        size_t result;
        char path[100];
        pFile = fopen(getFilePath(name, type, path), "r");
        fseek(pFile, 0, SEEK_END);
        lSize = ftell (pFile);
        rewind (pFile);
        
        // allocate memory to contain the whole file:
        buffer = (char*) malloc (sizeof(char)*lSize);
        if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}
        
        // copy the file into the buffer:
        result = fread (buffer,1,lSize,pFile);
        if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
        
        
        data_len = (size_t)lSize;
        if (buffer == NULL) {
            exit(-1);
            /* return 0; */
        }
        printf("filesize: %d\n", (int)data_len);
        
        {
            unsigned int flags = TINYOBJ_FLAG_TRIANGULATE;
            int ret = tinyobj_parse_obj(&attrib, &shapes, &num_shapes, &materials,
                                        &num_materials, buffer, data_len, flags);
            if (ret != TINYOBJ_SUCCESS) {
                return 0;
            }
            
            printf("# of shapes    = %d\n", (int)num_shapes);
            printf("# of materiasl = %d\n", (int)num_materials);
            
            /*
             {
             int i;
             for (i = 0; i < num_shapes; i++) {
             printf("shape[%d] name = %s\n", i, shapes[i].name);
             }
             }
             */
        }
        
        bmin[0] = bmin[1] = bmin[2] = FLT_MAX;
        bmax[0] = bmax[1] = bmax[2] = -FLT_MAX;
        
        {
            //DrawObject o;
            float* vb;
            /* std::vector<float> vb; //  */
            size_t face_offset = 0;
            size_t i;
            
            /* Assume triangulated face. */
            size_t num_triangles = attrib.num_face_num_verts;
            size_t stride = 9; /* 9 = pos(3float), normal(3float), color(3float) */
            
            vb = (float*)malloc(sizeof(float) * stride * num_triangles * 3);
            
            for (i = 0; i < attrib.num_face_num_verts; i++) {
                size_t f;
                assert(attrib.face_num_verts[i] % 3 == 0); /* assume all triangle faces. */
                for (f = 0; f < (size_t)attrib.face_num_verts[i] / 3; f++) {
                    size_t k;
                    float v[3][3];
                    float n[3][3];
                    float c[3];
                    float len2;
                    
                    tinyobj_vertex_index_t idx0 = attrib.faces[face_offset + 3 * f + 0];
                    tinyobj_vertex_index_t idx1 = attrib.faces[face_offset + 3 * f + 1];
                    tinyobj_vertex_index_t idx2 = attrib.faces[face_offset + 3 * f + 2];
                    
                    for (k = 0; k < 3; k++) {
                        int f0 = idx0.v_idx;
                        int f1 = idx1.v_idx;
                        int f2 = idx2.v_idx;
                        assert(f0 >= 0);
                        assert(f1 >= 0);
                        assert(f2 >= 0);
                        
                        v[0][k] = attrib.vertices[3 * (size_t)f0 + k];
                        v[1][k] = attrib.vertices[3 * (size_t)f1 + k];
                        v[2][k] = attrib.vertices[3 * (size_t)f2 + k];
                        bmin[k] = (v[0][k] < bmin[k]) ? v[0][k] : bmin[k];
                        bmin[k] = (v[1][k] < bmin[k]) ? v[1][k] : bmin[k];
                        bmin[k] = (v[2][k] < bmin[k]) ? v[2][k] : bmin[k];
                        bmax[k] = (v[0][k] > bmax[k]) ? v[0][k] : bmax[k];
                        bmax[k] = (v[1][k] > bmax[k]) ? v[1][k] : bmax[k];
                        bmax[k] = (v[2][k] > bmax[k]) ? v[2][k] : bmax[k];
                    }
                    
                    if (attrib.num_normals > 0) {
                        int f0 = idx0.vn_idx;
                        int f1 = idx1.vn_idx;
                        int f2 = idx2.vn_idx;
                        if (f0 >= 0 && f1 >= 0 && f2 >= 0) {
                            assert(f0 < (int)attrib.num_normals);
                            assert(f1 < (int)attrib.num_normals);
                            assert(f2 < (int)attrib.num_normals);
                            for (k = 0; k < 3; k++) {
                                n[0][k] = attrib.normals[3 * (size_t)f0 + k];
                                n[1][k] = attrib.normals[3 * (size_t)f1 + k];
                                n[2][k] = attrib.normals[3 * (size_t)f2 + k];
                            }
                        } else { /* normal index is not defined for this face */
                            /* compute geometric normal */
                            CalcNormal(n[0], v[0], v[1], v[2]);
                            n[1][0] = n[0][0];
                            n[1][1] = n[0][1];
                            n[1][2] = n[0][2];
                            n[2][0] = n[0][0];
                            n[2][1] = n[0][1];
                            n[2][2] = n[0][2];
                        }
                    } else {
                        /* compute geometric normal */
                        CalcNormal(n[0], v[0], v[1], v[2]);
                        n[1][0] = n[0][0];
                        n[1][1] = n[0][1];
                        n[1][2] = n[0][2];
                        n[2][0] = n[0][0];
                        n[2][1] = n[0][1];
                        n[2][2] = n[0][2];
                    }
                    
                    for (k = 0; k < 3; k++) {
                        vb[(3 * i + k) * stride + 0] = v[k][0];
                        vb[(3 * i + k) * stride + 1] = v[k][1];
                        vb[(3 * i + k) * stride + 2] = v[k][2];
                        vb[(3 * i + k) * stride + 3] = n[k][0];
                        vb[(3 * i + k) * stride + 4] = n[k][1];
                        vb[(3 * i + k) * stride + 5] = n[k][2];
                        
                        /* Use normal as color. */
                        c[0] = n[k][0];
                        c[1] = n[k][1];
                        c[2] = n[k][2];
                        len2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
                        if (len2 > 0.0f) {
                            float len = (float)sqrt((double)len2);
                            
                            c[0] /= len;
                            c[1] /= len;
                            c[2] /= len;
                        }
                        
                        vb[(3 * i + k) * stride + 6] = (c[0] * 0.5f + 0.5f);
                        vb[(3 * i + k) * stride + 7] = (c[1] * 0.5f + 0.5f);
                        vb[(3 * i + k) * stride + 8] = (c[2] * 0.5f + 0.5f);
                    }
                }
                face_offset += (size_t)attrib.face_num_verts[i];
            }
            
            //o.vb = 0;
            //o.numTriangles = 0;
            //        if (num_triangles > 0) {
            //            glGenBuffers(1, &o.vb);
            //            glBindBuffer(GL_ARRAY_BUFFER, o.vb);
            //            glBufferData(GL_ARRAY_BUFFER, num_triangles * 3 * stride * sizeof(float),
            //                         vb, GL_STATIC_DRAW);
            //            o.numTriangles = (int)num_triangles;
            //        }
            
            free(vb);
            
            //gDrawObject = o;
        }
        
        printf("bmin = %f, %f, %f\n", (double)bmin[0], (double)bmin[1],
               (double)bmin[2]);
        printf("bmax = %f, %f, %f\n", (double)bmax[0], (double)bmax[1],
               (double)bmax[2]);
        
        tinyobj_attrib_free(&attrib);
        tinyobj_shapes_free(shapes, num_shapes);
        tinyobj_materials_free(materials, num_materials);
        
        // terminate
        fclose (pFile);
        free (buffer);
    }
    return 0;
}

#define PNG_BYTES_TO_CHECK 4
int load_png_image( const char *name, const char *type, unsigned int **bits, unsigned int *width, unsigned int *height)
{
    FILE *fp;
    png_structp png_ptr;
    png_infop info_ptr;
    png_bytep* row_pointers;
    char buf[PNG_BYTES_TO_CHECK];
    int w, h, x, y, temp, color_type;
    
    char path[100];
    fp = fopen( getFilePath(name, type, path), "rb" );
    if( fp == NULL ) {
        return 1; /* 返回值 */
    }
    
    png_ptr = png_create_read_struct( PNG_LIBPNG_VER_STRING, 0, 0, 0 );
    info_ptr = png_create_info_struct( png_ptr );
    
    setjmp( png_jmpbuf(png_ptr) );
    /* 读取PNG_BYTES_TO_CHECK个字节的数据 */
    temp = (int)fread( buf, 1, PNG_BYTES_TO_CHECK, fp );
    /* 若读到的数据并没有PNG_BYTES_TO_CHECK个字节 */
    if( temp < PNG_BYTES_TO_CHECK ) {
        fclose(fp);
        png_destroy_read_struct( &png_ptr, &info_ptr, 0);
        return 2;/* 返回值 */
    }
    /* 检测数据是否为PNG的签名 */
    temp = png_sig_cmp( (png_bytep)buf, (png_size_t)0, PNG_BYTES_TO_CHECK );
    /* 如果不是PNG的签名，则说明该文件不是PNG文件 */
    if( temp != 0 ) {
        fclose(fp);
        png_destroy_read_struct( &png_ptr, &info_ptr, 0);
        return 3;/* 返回值 */
    }
    
    /* 复位文件指针 */
    rewind( fp );
    /* 开始读文件 */
    png_init_io( png_ptr, fp );
    /* 读取PNG图片信息 */
    png_read_png( png_ptr, info_ptr, PNG_TRANSFORM_EXPAND, 0 );
    /* 获取图像的色彩类型 */
    color_type = png_get_color_type( png_ptr, info_ptr );
    /* 获取图像的宽高 */
    w = png_get_image_width( png_ptr, info_ptr );
    h = png_get_image_height( png_ptr, info_ptr );
    
    *bits = (unsigned int*)malloc(sizeof(unsigned int) * w * h);
    
    /* 获取图像的所有行像素数据，row_pointers里边就是rgba数据 */
    row_pointers = png_get_rows( png_ptr, info_ptr );
    /* 根据不同的色彩类型进行相应处理 */
    switch( color_type ) {
        case PNG_COLOR_TYPE_RGB_ALPHA:
            for( y=0; y<h; ++y ) {
                for( x=0; x < w; ++x ) {
                    (*bits)[y*w+x] = 0;
                    /* 以下是RGBA数据，需要自己补充代码，保存RGBA数据 */
                    (*bits)[y*w+x] |= row_pointers[y][4*x+0] << 16; // red
                    (*bits)[y*w+x] |= row_pointers[y][4*x+1] << 8; // green
                    (*bits)[y*w+x] |= row_pointers[y][4*x+2]; // blue
                    (*bits)[y*w+x] |= row_pointers[y][4*x+3] << 24; // alpha
                }
            }
            break;
            
        case PNG_COLOR_TYPE_RGB:
            for( y=0; y<h; ++y ) {
                for( x=0; x<w; ++x ) {
                    (*bits)[y*w+x] = 0xff000000;
                    (*bits)[y*w+x] |= row_pointers[y][3*x+0] << 16; // red
                    (*bits)[y*w+x] |= row_pointers[y][3*x+1] << 8; // green
                    (*bits)[y*w+x] |= row_pointers[y][3*x+2]; // blue
                }
            }
            break;
            /* 其它色彩类型的图像就不读了 */
        default:
            fclose(fp);
            png_destroy_read_struct( &png_ptr, &info_ptr, 0);
            return 4/* 返回值 */;
    }
    png_destroy_read_struct( &png_ptr, &info_ptr, 0);
    
    *width = w;
    *height = h;
    
    return 0;
}
