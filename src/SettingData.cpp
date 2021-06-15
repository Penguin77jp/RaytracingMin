#include "SettingData.h"
#include <fstream>

namespace png {
	LoadData::LoadData(std::string jsonName) {
		std::ifstream jsonStream(jsonName);
		nlohmann::json json;
		jsonStream >> json;

		from_json(json, data);
	}

	void LoadData::SaveSampleJson(std::string fileName) {
		SettingData tmp;
		tmp.width = 960;
		tmp.height = 540;
		tmp.sample = 5;
		tmp.camera.origin = vec3(0, 0, 0);
		tmp.camera.target = vec3(0, 0, 1);
		tmp.camera.upVec = vec3(0, 1, 0);
		tmp.camera.fov = 60;

		tmp.object.push_back({
			vec3(0,0,10)
			, 1.0
			, Material{
				vec3(1.0,0.2,0.2)
				,0.8
				,vec3(0,0,0)
			}
			});
		//light
		tmp.object.push_back({
			vec3(0,5,0)
			, 1.0
			, Material{
				vec3(1.0,1.0,1.0)
				,1.0
				,vec3(1,1,1)
			}
			});
		//wall
		{
			float kd = 0.8;
			float size = 99999;
			float space = 10;
			tmp.object.push_back({
				vec3(size+0.5*space,0,0)
				, size
				, Material{
					vec3(1.0,0.5,0.5)
					,kd
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(-(size + 0.5 * space),0,0)
				, size
				, Material{
					vec3(0.5,1.0,0.5)
					,kd
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,size + 0.5 * space,0)
				, size
				, Material{
					vec3(0.5,0.5,1.0)
					,kd
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,-(size + 0.5 * space),0)
				, size
				, Material{
					vec3(0.5,1.0,1.0)
					,kd
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,0,size + 0.5 * space)
				, size
				, Material{
					vec3(1.0,0.5,1.0)
					,kd
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,0,-(size + 0.5 * space))
				, size
				, Material{
					vec3(1.0,1.0,0.5)
					,kd
					,vec3()
				}
				});
		}

		nlohmann::json tmpJson;
		to_json(tmpJson, tmp);
		std::ofstream ostream(fileName);
		ostream << tmpJson;
	}

	void LoadData::to_json(nlohmann::json& json, const SettingData& data) {
		json = {
			{"00 width", data.width}
			,{"00 height", data.height}
			,{"00 sample",data.sample}
			,{"01 camera",{
				{"origin",{data.camera.origin.x,data.camera.origin.y,data.camera.origin.z}}
				,{"target",{data.camera.target.x,data.camera.target.y,data.camera.target.z}}
				,{"upVec",{data.camera.upVec.x,data.camera.upVec.y,data.camera.upVec.z}}
				,{"fov",data.camera.fov}
			}}
		};
		for (int i = 0; i < data.object.size(); ++i) {
			auto& obj = data.object[i];
			json["02 scene"]["00 object"][i]["00 position"] = { obj.position.x,obj.position.y,obj.position.z };
			json["02 scene"]["00 object"][i]["01 size"] = obj.size;
			json["02 scene"]["00 object"][i]["02 material"]["color"] = { obj.material.color.x ,obj.material.color.y,obj.material.color.z };
			json["02 scene"]["00 object"][i]["02 material"]["kd"] = obj.material.kd;
			json["02 scene"]["00 object"][i]["02 material"]["emission"] = { obj.material.emission.x,obj.material.emission.y,obj.material.emission.z };
		}
	}
	void LoadData::from_json(const nlohmann::json& json, SettingData& data) {
		for (auto& it : json.items()) {
			if (it.key() == "00 width") data.width = it.value();
			else if (it.key() == "00 height") data.height = it.value();
			else if (it.key() == "00 sample") data.sample = it.value();
			else if (it.key() == "01 camera") {
				for (auto& it_cam : it.value().items()) {
					if (it_cam.key() == "origin") {
						data.camera.origin.x = it_cam.value()[0];
						data.camera.origin.y = it_cam.value()[1];
						data.camera.origin.z = it_cam.value()[2];
					}
					else if (it_cam.key() == "target") {
						data.camera.target.x = it_cam.value()[0];
						data.camera.target.y = it_cam.value()[1];
						data.camera.target.z = it_cam.value()[2];
					}
					else if (it_cam.key() == "upVec") {
						data.camera.upVec.x = it_cam.value()[0];
						data.camera.upVec.y = it_cam.value()[1];
						data.camera.upVec.z = it_cam.value()[2];
					}
					else if (it_cam.key() == "fov") {
						data.camera.fov = it_cam.value();
					}
				}
			}
			else if (it.key() == "02 scene") {
				for (auto& it_scene : it.value().items()) {
					if (it_scene.key() == "00 object") {
						for (auto& it_object : it_scene.value().items()) {
							Object tmp;
							tmp.position.x = it_object.value()["00 position"][0];
							tmp.position.y = it_object.value()["00 position"][1];
							tmp.position.z = it_object.value()["00 position"][2];
							tmp.size = it_object.value()["01 size"];
							tmp.material.color.x = it_object.value()["02 material"]["color"][0];
							tmp.material.color.y = it_object.value()["02 material"]["color"][1];
							tmp.material.color.z = it_object.value()["02 material"]["color"][2];
							tmp.material.kd = it_object.value()["02 material"]["kd"];
							tmp.material.emission.x = it_object.value()["02 material"]["emission"][0];
							tmp.material.emission.y = it_object.value()["02 material"]["emission"][1];
							tmp.material.emission.z = it_object.value()["02 material"]["emission"][2];
							data.object.push_back(tmp);
						}
					}
				}
			}
		}

	}

}