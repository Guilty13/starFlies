#include <memory>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <evhttp.h>

#include "stars.h"

using namespace std;

string json(const vector<float>& xs) {
    string res = "[";
    for (size_t i = 0; i < xs.size(); i++) {
        if (i > 0) {
            res += ",";
        }
        res += to_string(xs[i]);
    }
    res += "]";
    return res;
}

//FIXME: macro
void VERIFY(bool p) {
    if (!p) {
        //FIXME: std::exception
        throw string("VEFRIFY");
    }
}

int main()
{
    if (!event_init())
    {
        std::cerr << "Failed to init libevent." << std::endl;
        return -1;
    }
    char const SrvAddress[] = "*"; //"127.0.0.1";
    std::uint16_t SrvPort = 28080;
    std::unique_ptr<evhttp, decltype(&evhttp_free)> Server(evhttp_start(SrvAddress, SrvPort), &evhttp_free);
    if (!Server)
    {
        std::cerr << "Failed to init http server." << std::endl;
        return -1;
    }
    auto OnReq = [] (evhttp_request *req, void *)
    {
        try {
            auto uri = evhttp_request_get_evhttp_uri(req);
            string path = evhttp_uri_get_path(uri);

            auto *OutBuf = evhttp_request_get_output_buffer(req);

            if (path == "/") {

                ifstream ifs("flightsNew.html"); //"main.html");
                string str((istreambuf_iterator<char>(ifs)), istreambuf_iterator<char>());
            
                evbuffer_add_printf(OutBuf, "%s", str.c_str());
            } else if (path == "/compute.jsonp") {

                evkeyvalq params;
                evhttp_parse_query_str(evhttp_uri_get_query(uri), &params);

                VERIFY(evhttp_find_header(&params, "x0"));
                VERIFY(evhttp_find_header(&params, "y0"));
                VERIFY(evhttp_find_header(&params, "vx0"));
                VERIFY(evhttp_find_header(&params, "vy0"));

                float x0 = atof(evhttp_find_header(&params, "x0"));
                float y0 = atof(evhttp_find_header(&params, "y0"));
                float vx0 = atof(evhttp_find_header(&params, "vx0"));
                float vy0 = atof(evhttp_find_header(&params, "vy0"));
                cerr << "X=" << x0 << " Y=" << y0 << " VX=" << vx0 << " VY=" << vy0 << endl;

                auto ctx = evhttp_find_header(&params, "ctx");
                auto route = Compute({x0,y0,vx0,vy0});

                string callback = evhttp_find_header(&params, "callback");
                string out;
                out += callback + "({" +
                    "'route': " + json(route) +
                    (ctx ? string(", 'ctx': ") + ctx : "") +
                "})\n";
                evbuffer_add_printf(OutBuf, out.c_str());
            }
            evhttp_send_reply(req, HTTP_OK, "", OutBuf);
        } catch (...) {
            evhttp_send_reply(req, 500, "WTF", nullptr);
        }
    };
    evhttp_set_gencb(Server.get(), OnReq, nullptr);
    if (event_dispatch() == -1)
    {
        std::cerr << "Failed to run messahe loop." << std::endl;
        return -1;
    }
    return 0;
}
