#include "SettingData.h"
#include <fstream>

namespace png {
	LoadData::LoadData(std::string jsonName) {
		std::ifstream jsonStream(jsonName);
		nlohmann::json json;
		jsonStream >> json;

		// innitialize strcut
		{
			data.renderType = 0;
			data.width = data.height = 100;
			data.samples = data.spectrumSamples = 10;
			data.superSamples = 1;
			data.cameraOrigin = vec3();
			data.cameraTarget = vec3(0, 0, 1);
			data.fov = 1.0;
		}

		//debug
		//data.object.push_back(new BoxObject(vec3(), vec3(1,1,1), new RefractionMaterial(TransparentMaterialType::BK7, new TextureSolid(vec3(1, 1, 1)), new TextureSolid(vec3()))));
		data.object.push_back(new WavePlane(vec3(0, 0, 0), nullptr));
		data.object.push_back(new SphereObject(vec3(0, 0.4, 1.5), 0.2, new DiffuseMaterial(new TextureSolid(vec3(0)), new TextureSolid(vec3(1, 1, 1)))));
		data.object.push_back(new SphereObject(vec3(0.5,0.4,1.5), 0.3, new DiffuseMaterial(new TextureSolid(vec3(0)), new TextureSolid(vec3(0,0,1)))));
		data.object.push_back(new SphereObject(vec3(0, 0.7, 0.0), 0.3, new DiffuseMaterial(new TextureSolid(vec3(0)), new TextureSolid(vec3(1,0,0)))));

		// glass
		//data.object.push_back(new SphereObject(vec3(0,0,0), 3.0, new RefractionMaterial(TransparentMaterialType::BK7, new TextureSolid(vec3(1,1,1)), new TextureSolid(vec3()))));

		// emistor image
		//data.object.push_back(new PlaneObject(vec3(0, 0, 9), vec3(0, 10, 0), vec3(10, 0, 0), new DiffuseMaterial(new TextureSolid(vec3()), new TextureImage("penguin.jpg", 2.0f))));

		from_json(json, data);
	}

	void LoadData::SaveSampleJson(std::string fileName) {
		SettingData tmp;
		tmp.renderType = 0;
		tmp.width = 960;
		tmp.height = 540;
		tmp.samples = 10;
		tmp.superSamples = 4;
		tmp.cameraOrigin = vec3(0, 0, 0);
		tmp.cameraTarget = vec3(0, 0, 1);
		tmp.fov = 60;

		/*
		tmp.object.push_back({
			vec3(0,0,10)
			, 1.0
			, Material{
				vec3(1.0,0.2,0.2)
				,vec3(0,0,0)
			}
			});

		//light
		tmp.object.push_back({
			vec3(0,11,10)
			, 10.0
			, Material{
				vec3(1.0,1.0,1.0)
				,vec3(1,1,1)
			}
			});

		//wall
		{
			float size = 99999;
			float space = 10;
			tmp.object.push_back({
				vec3(size+0.5*space,0,0)
				, size
				, Material{
					vec3(1.0,0.5,0.5)
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(-(size + 0.5 * space),0,0)
				, size
				, Material{
					vec3(0.5,1.0,0.5)
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,size + 0.5 * space,0)
				, size
				, Material{
					vec3(0.5,0.5,1.0)
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,-(size + 0.5 * space),0)
				, size
				, Material{
					vec3(0.5,1.0,1.0)
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,0,size + 0.5 * space)
				, size
				, Material{
					vec3(1.0,0.5,1.0)
					,vec3()
				}
				});
			tmp.object.push_back({
				vec3(0,0,-(size + 0.5 * space))
				, size
				, Material{
					vec3(1.0,1.0,0.5)
					,vec3()
				}
				});
		}
		*/

		nlohmann::json tmpJson;
		to_json(tmpJson, tmp);
		std::ofstream ostream(fileName);
		ostream << tmpJson;
	}

	void LoadData::to_json(nlohmann::json& json, const SettingData& data) {
		json = {
			{"00 width", data.width}
			,{"00 height", data.height}
			,{"00 samples",data.samples}
			,{"00 superSamples",data.superSamples}
			,{"01 camera",{
				{"origin",{data.cameraOrigin.x,data.cameraOrigin.y,data.cameraOrigin.z}}
				,{"target",{data.cameraTarget.x,data.cameraTarget.y,data.cameraTarget.z}}
				,{"fov",data.fov}
			}}
		};
		for (int i = 0; i < data.object.size(); ++i) {
			auto& obj = data.object[i];
			/*
			json["02 scene"]["00 object"][i]["00 position"] = { obj.position.x,obj.position.y,obj.position.z };
			json["02 scene"]["00 object"][i]["01 size"] = obj.size;
			json["02 scene"]["00 object"][i]["02 material"]["color"] = { obj.material.color.x ,obj.material.color.y,obj.material.color.z };
			json["02 scene"]["00 object"][i]["02 material"]["emission"] = { obj.material.emission.x,obj.material.emission.y,obj.material.emission.z };
			*/
		}
	}
	void LoadData::from_json(const nlohmann::json& json, SettingData& data) {
		for (auto& it : json.items()) {
			if (it.key() == "00 renderType") data.renderType = it.value();
			else if (it.key() == "00 width") data.width = it.value();
			else if (it.key() == "00 height") data.height = it.value();
			else if (it.key() == "00 samples") data.samples = it.value();
			else if (it.key() == "00 superSamples") data.superSamples = it.value();
			else if (it.key() == "00 spectrumSamples") data.spectrumSamples = it.value();
			else if (it.key() == "01 camera") {
				for (auto& it_cam : it.value().items()) {
					if (it_cam.key() == "origin") {
						data.cameraOrigin.x = it_cam.value()[0];
						data.cameraOrigin.y = it_cam.value()[1];
						data.cameraOrigin.z = it_cam.value()[2];
					}
					else if (it_cam.key() == "target") {
						data.cameraTarget.x = it_cam.value()[0];
						data.cameraTarget.y = it_cam.value()[1];
						data.cameraTarget.z = it_cam.value()[2];
					}
					else if (it_cam.key() == "fov") {
						data.fov = it_cam.value();
					}
				}
			}
			else if (it.key() == "02 scene") {
				for (auto& it_scene : it.value().items()) {
					if (it_scene.key() == "00 object") {
						for (auto& it_object : it_scene.value().items()) {
							SceneObject* tmp;
							int objectType = it_object.value()["00 objectType"];
							// sphere
							if (objectType == 0) {
								vec3 posi;
								posi.x = it_object.value()["01 position"][0];
								posi.y = it_object.value()["01 position"][1];
								posi.z = it_object.value()["01 position"][2];
								float size = it_object.value()["02 size"];

								//material
								int materialType = it_object.value()["03 material"]["type"];
								Material* mat = nullptr;
								if (materialType == 0) {
									vec3 color;
									color.x = it_object.value()["03 material"]["color"][0];
									color.y = it_object.value()["03 material"]["color"][1];
									color.z = it_object.value()["03 material"]["color"][2];
									vec3 emission;
									emission.x = it_object.value()["03 material"]["emission"][0];
									emission.y = it_object.value()["03 material"]["emission"][1];
									emission.z = it_object.value()["03 material"]["emission"][2];
									mat = new DiffuseMaterial(new TextureSolid(color), new TextureSolid(emission));
								}
								else if (materialType == 1) {
									vec3 color;
									color.x = it_object.value()["03 material"]["color"][0];
									color.y = it_object.value()["03 material"]["color"][1];
									color.z = it_object.value()["03 material"]["color"][2];
									vec3 emission;
									emission.x = it_object.value()["03 material"]["emission"][0];
									emission.y = it_object.value()["03 material"]["emission"][1];
									emission.z = it_object.value()["03 material"]["emission"][2];
									mat = new RefractionMaterial(TransparentMaterialType::HighVariance, new TextureSolid(color), new TextureSolid(emission));
								}

								tmp = new SphereObject(posi, size, mat);
							}

							// wave plane
							else if (objectType == 1) {

								tmp = new WavePlane(vec3(0,0,0), nullptr);
							}

							data.object.push_back(tmp);
						}
					}
				}
			}
		}

	}

}