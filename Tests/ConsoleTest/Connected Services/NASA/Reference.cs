﻿//------------------------------------------------------------------------------
// <auto-generated>
//     Этот код создан программой.
//
//     Изменения в этом файле могут привести к неправильной работе и будут потеряны в случае
//     повторного создания кода.
// </auto-generated>
//------------------------------------------------------------------------------

namespace NASA
{
    
    
    /// <remarks/>
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.Xml.Serialization.XmlTypeAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/")]
    public partial class HelioExternalException
    {
        
        private string messageField;
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=0)]
        public string message
        {
            get
            {
                return this.messageField;
            }
            set
            {
                this.messageField = value;
            }
        }
    }
    
    /// <remarks/>
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.Xml.Serialization.XmlTypeAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/")]
    public partial class objectDescription
    {
        
        private System.DateTime endDateField;
        
        private bool endDateFieldSpecified;
        
        private string idField;
        
        private string nameField;
        
        private System.DateTime startDateField;
        
        private bool startDateFieldSpecified;
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=0)]
        public System.DateTime endDate
        {
            get
            {
                return this.endDateField;
            }
            set
            {
                this.endDateField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlIgnoreAttribute()]
        public bool endDateSpecified
        {
            get
            {
                return this.endDateFieldSpecified;
            }
            set
            {
                this.endDateFieldSpecified = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=1)]
        public string id
        {
            get
            {
                return this.idField;
            }
            set
            {
                this.idField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=2)]
        public string name
        {
            get
            {
                return this.nameField;
            }
            set
            {
                this.nameField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=3)]
        public System.DateTime startDate
        {
            get
            {
                return this.startDateField;
            }
            set
            {
                this.startDateField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlIgnoreAttribute()]
        public bool startDateSpecified
        {
            get
            {
                return this.startDateFieldSpecified;
            }
            set
            {
                this.startDateFieldSpecified = value;
            }
        }
    }
    
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.ServiceModel.ServiceContractAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/", ConfigurationName="NASA.HeliocentricTrajectoriesInterface")]
    public interface HeliocentricTrajectoriesInterface
    {
        
        [System.ServiceModel.OperationContractAttribute(Action="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getAllObjectsRe" +
            "quest", ReplyAction="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getAllObjectsRe" +
            "sponse")]
        [System.ServiceModel.FaultContractAttribute(typeof(NASA.HelioExternalException), Action="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getAllObjects/F" +
            "ault/HelioExternalException", Name="HelioExternalException")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<NASA.getAllObjectsResponse> getAllObjectsAsync(NASA.getAllObjectsRequest request);
        
        [System.ServiceModel.OperationContractAttribute(Action="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getTrajectories" +
            "Request", ReplyAction="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getTrajectories" +
            "Response")]
        [System.ServiceModel.FaultContractAttribute(typeof(NASA.HelioExternalException), Action="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getTrajectories" +
            "/Fault/HelioExternalException", Name="HelioExternalException")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<NASA.getTrajectoriesResponse> getTrajectoriesAsync(NASA.getTrajectoriesRequest request);
        
        [System.ServiceModel.OperationContractAttribute(Action="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getPrivacyAndIm" +
            "portantNoticesRequest", ReplyAction="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getPrivacyAndIm" +
            "portantNoticesResponse")]
        [System.ServiceModel.FaultContractAttribute(typeof(NASA.HelioExternalException), Action="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getPrivacyAndIm" +
            "portantNotices/Fault/HelioExternalException", Name="HelioExternalException")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<NASA.getPrivacyAndImportantNoticesResponse> getPrivacyAndImportantNoticesAsync(NASA.getPrivacyAndImportantNoticesRequest request);
        
        [System.ServiceModel.OperationContractAttribute(Action="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getAcknowledgem" +
            "entsRequest", ReplyAction="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getAcknowledgem" +
            "entsResponse")]
        [System.ServiceModel.FaultContractAttribute(typeof(NASA.HelioExternalException), Action="http://helio.spdf.gsfc.nasa.gov/HeliocentricTrajectoriesInterface/getAcknowledgem" +
            "ents/Fault/HelioExternalException", Name="HelioExternalException")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<NASA.getAcknowledgementsResponse> getAcknowledgementsAsync(NASA.getAcknowledgementsRequest request);
    }
    
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
    [System.ServiceModel.MessageContractAttribute(WrapperName="getAllObjects", WrapperNamespace="http://helio.spdf.gsfc.nasa.gov/", IsWrapped=true)]
    public partial class getAllObjectsRequest
    {
        
        public getAllObjectsRequest()
        {
        }
    }
    
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
    [System.ServiceModel.MessageContractAttribute(WrapperName="getAllObjectsResponse", WrapperNamespace="http://helio.spdf.gsfc.nasa.gov/", IsWrapped=true)]
    public partial class getAllObjectsResponse
    {
        
        [System.ServiceModel.MessageBodyMemberAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/", Order=0)]
        [System.Xml.Serialization.XmlElementAttribute("return", Form=System.Xml.Schema.XmlSchemaForm.Unqualified)]
        public NASA.objectDescription[] @return;
        
        public getAllObjectsResponse()
        {
        }
        
        public getAllObjectsResponse(NASA.objectDescription[] @return)
        {
            this.@return = @return;
        }
    }
    
    /// <remarks/>
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.Xml.Serialization.XmlTypeAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/")]
    public partial class request
    {
        
        private System.DateTime beginDateField;
        
        private bool beginDateFieldSpecified;
        
        private coordinateSystem coordinateSystemField;
        
        private bool coordinateSystemFieldSpecified;
        
        private System.DateTime endDateField;
        
        private bool endDateFieldSpecified;
        
        private string[] objectsField;
        
        private int resolutionField;
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=0)]
        public System.DateTime beginDate
        {
            get
            {
                return this.beginDateField;
            }
            set
            {
                this.beginDateField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlIgnoreAttribute()]
        public bool beginDateSpecified
        {
            get
            {
                return this.beginDateFieldSpecified;
            }
            set
            {
                this.beginDateFieldSpecified = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=1)]
        public coordinateSystem coordinateSystem
        {
            get
            {
                return this.coordinateSystemField;
            }
            set
            {
                this.coordinateSystemField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlIgnoreAttribute()]
        public bool coordinateSystemSpecified
        {
            get
            {
                return this.coordinateSystemFieldSpecified;
            }
            set
            {
                this.coordinateSystemFieldSpecified = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=2)]
        public System.DateTime endDate
        {
            get
            {
                return this.endDateField;
            }
            set
            {
                this.endDateField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlIgnoreAttribute()]
        public bool endDateSpecified
        {
            get
            {
                return this.endDateFieldSpecified;
            }
            set
            {
                this.endDateFieldSpecified = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute("objects", Form=System.Xml.Schema.XmlSchemaForm.Unqualified, IsNullable=true, Order=3)]
        public string[] objects
        {
            get
            {
                return this.objectsField;
            }
            set
            {
                this.objectsField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=4)]
        public int resolution
        {
            get
            {
                return this.resolutionField;
            }
            set
            {
                this.resolutionField = value;
            }
        }
    }
    
    /// <remarks/>
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.Xml.Serialization.XmlTypeAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/")]
    public enum coordinateSystem
    {
        
        /// <remarks/>
        SE,
        
        /// <remarks/>
        HG,
        
        /// <remarks/>
        HGI,
    }
    
    /// <remarks/>
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.Xml.Serialization.XmlTypeAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/")]
    public partial class dataResult : result
    {
        
        private trajectory[] trajectoriesField;
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute("trajectories", Form=System.Xml.Schema.XmlSchemaForm.Unqualified, IsNullable=true, Order=0)]
        public trajectory[] trajectories
        {
            get
            {
                return this.trajectoriesField;
            }
            set
            {
                this.trajectoriesField = value;
            }
        }
    }
    
    /// <remarks/>
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.Xml.Serialization.XmlTypeAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/")]
    public partial class trajectory
    {
        
        private coordinateSystem coordinateSystemField;
        
        private bool coordinateSystemFieldSpecified;
        
        private string idField;
        
        private System.Nullable<double>[] latitudeField;
        
        private System.Nullable<double>[] longitudeField;
        
        private System.Nullable<double>[] radiusField;
        
        private System.Nullable<System.DateTime>[] timeField;
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=0)]
        public coordinateSystem coordinateSystem
        {
            get
            {
                return this.coordinateSystemField;
            }
            set
            {
                this.coordinateSystemField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlIgnoreAttribute()]
        public bool coordinateSystemSpecified
        {
            get
            {
                return this.coordinateSystemFieldSpecified;
            }
            set
            {
                this.coordinateSystemFieldSpecified = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=1)]
        public string id
        {
            get
            {
                return this.idField;
            }
            set
            {
                this.idField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute("latitude", Form=System.Xml.Schema.XmlSchemaForm.Unqualified, IsNullable=true, Order=2)]
        public System.Nullable<double>[] latitude
        {
            get
            {
                return this.latitudeField;
            }
            set
            {
                this.latitudeField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute("longitude", Form=System.Xml.Schema.XmlSchemaForm.Unqualified, IsNullable=true, Order=3)]
        public System.Nullable<double>[] longitude
        {
            get
            {
                return this.longitudeField;
            }
            set
            {
                this.longitudeField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute("radius", Form=System.Xml.Schema.XmlSchemaForm.Unqualified, IsNullable=true, Order=4)]
        public System.Nullable<double>[] radius
        {
            get
            {
                return this.radiusField;
            }
            set
            {
                this.radiusField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute("time", Form=System.Xml.Schema.XmlSchemaForm.Unqualified, IsNullable=true, Order=5)]
        public System.Nullable<System.DateTime>[] time
        {
            get
            {
                return this.timeField;
            }
            set
            {
                this.timeField = value;
            }
        }
    }
    
    /// <remarks/>
    [System.Xml.Serialization.XmlIncludeAttribute(typeof(dataResult))]
    [System.Xml.Serialization.XmlIncludeAttribute(typeof(fileResult))]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.Xml.Serialization.XmlTypeAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/")]
    public partial class result
    {
        
        private resultStatusCode statusCodeField;
        
        private bool statusCodeFieldSpecified;
        
        private resultStatusSubCode statusSubCodeField;
        
        private bool statusSubCodeFieldSpecified;
        
        private string[] statusTextField;
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=0)]
        public resultStatusCode statusCode
        {
            get
            {
                return this.statusCodeField;
            }
            set
            {
                this.statusCodeField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlIgnoreAttribute()]
        public bool statusCodeSpecified
        {
            get
            {
                return this.statusCodeFieldSpecified;
            }
            set
            {
                this.statusCodeFieldSpecified = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified, Order=1)]
        public resultStatusSubCode statusSubCode
        {
            get
            {
                return this.statusSubCodeField;
            }
            set
            {
                this.statusSubCodeField = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlIgnoreAttribute()]
        public bool statusSubCodeSpecified
        {
            get
            {
                return this.statusSubCodeFieldSpecified;
            }
            set
            {
                this.statusSubCodeFieldSpecified = value;
            }
        }
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute("statusText", Form=System.Xml.Schema.XmlSchemaForm.Unqualified, IsNullable=true, Order=2)]
        public string[] statusText
        {
            get
            {
                return this.statusTextField;
            }
            set
            {
                this.statusTextField = value;
            }
        }
    }
    
    /// <remarks/>
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.Xml.Serialization.XmlTypeAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/")]
    public enum resultStatusCode
    {
        
        /// <remarks/>
        SUCCESS,
        
        /// <remarks/>
        CONDITIONAL_SUCCESS,
        
        /// <remarks/>
        ERROR,
    }
    
    /// <remarks/>
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.Xml.Serialization.XmlTypeAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/")]
    public enum resultStatusSubCode
    {
        
        /// <remarks/>
        SUCCESS,
        
        /// <remarks/>
        INVALID_OBJ_ID,
        
        /// <remarks/>
        INVALID_TIME_RANGE,
        
        /// <remarks/>
        UPSTREAM_SERVER_ERROR,
        
        /// <remarks/>
        SERVER_ERROR,
    }
    
    /// <remarks/>
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.Xml.Serialization.XmlTypeAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/")]
    public partial class fileResult : result
    {
        
        private string[] urlsField;
        
        /// <remarks/>
        [System.Xml.Serialization.XmlElementAttribute("urls", Form=System.Xml.Schema.XmlSchemaForm.Unqualified, IsNullable=true, Order=0)]
        public string[] urls
        {
            get
            {
                return this.urlsField;
            }
            set
            {
                this.urlsField = value;
            }
        }
    }
    
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
    [System.ServiceModel.MessageContractAttribute(WrapperName="getTrajectories", WrapperNamespace="http://helio.spdf.gsfc.nasa.gov/", IsWrapped=true)]
    public partial class getTrajectoriesRequest
    {
        
        [System.ServiceModel.MessageBodyMemberAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/", Order=0)]
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified)]
        public NASA.request arg0;
        
        public getTrajectoriesRequest()
        {
        }
        
        public getTrajectoriesRequest(NASA.request arg0)
        {
            this.arg0 = arg0;
        }
    }
    
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
    [System.ServiceModel.MessageContractAttribute(WrapperName="getTrajectoriesResponse", WrapperNamespace="http://helio.spdf.gsfc.nasa.gov/", IsWrapped=true)]
    public partial class getTrajectoriesResponse
    {
        
        [System.ServiceModel.MessageBodyMemberAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/", Order=0)]
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified)]
        public NASA.dataResult @return;
        
        public getTrajectoriesResponse()
        {
        }
        
        public getTrajectoriesResponse(NASA.dataResult @return)
        {
            this.@return = @return;
        }
    }
    
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
    [System.ServiceModel.MessageContractAttribute(WrapperName="getPrivacyAndImportantNotices", WrapperNamespace="http://helio.spdf.gsfc.nasa.gov/", IsWrapped=true)]
    public partial class getPrivacyAndImportantNoticesRequest
    {
        
        public getPrivacyAndImportantNoticesRequest()
        {
        }
    }
    
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
    [System.ServiceModel.MessageContractAttribute(WrapperName="getPrivacyAndImportantNoticesResponse", WrapperNamespace="http://helio.spdf.gsfc.nasa.gov/", IsWrapped=true)]
    public partial class getPrivacyAndImportantNoticesResponse
    {
        
        [System.ServiceModel.MessageBodyMemberAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/", Order=0)]
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified)]
        public NASA.fileResult @return;
        
        public getPrivacyAndImportantNoticesResponse()
        {
        }
        
        public getPrivacyAndImportantNoticesResponse(NASA.fileResult @return)
        {
            this.@return = @return;
        }
    }
    
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
    [System.ServiceModel.MessageContractAttribute(WrapperName="getAcknowledgements", WrapperNamespace="http://helio.spdf.gsfc.nasa.gov/", IsWrapped=true)]
    public partial class getAcknowledgementsRequest
    {
        
        public getAcknowledgementsRequest()
        {
        }
    }
    
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
    [System.ServiceModel.MessageContractAttribute(WrapperName="getAcknowledgementsResponse", WrapperNamespace="http://helio.spdf.gsfc.nasa.gov/", IsWrapped=true)]
    public partial class getAcknowledgementsResponse
    {
        
        [System.ServiceModel.MessageBodyMemberAttribute(Namespace="http://helio.spdf.gsfc.nasa.gov/", Order=0)]
        [System.Xml.Serialization.XmlElementAttribute(Form=System.Xml.Schema.XmlSchemaForm.Unqualified)]
        public NASA.fileResult @return;
        
        public getAcknowledgementsResponse()
        {
        }
        
        public getAcknowledgementsResponse(NASA.fileResult @return)
        {
            this.@return = @return;
        }
    }
    
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    public interface HeliocentricTrajectoriesInterfaceChannel : NASA.HeliocentricTrajectoriesInterface, System.ServiceModel.IClientChannel
    {
    }
    
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    public partial class HeliocentricTrajectoriesInterfaceClient : System.ServiceModel.ClientBase<NASA.HeliocentricTrajectoriesInterface>, NASA.HeliocentricTrajectoriesInterface
    {
        
        /// <summary>
        /// Реализуйте этот разделяемый метод для настройки конечной точки службы.
        /// </summary>
        /// <param name="serviceEndpoint">Настраиваемая конечная точка</param>
        /// <param name="clientCredentials">Учетные данные клиента.</param>
        static partial void ConfigureEndpoint(System.ServiceModel.Description.ServiceEndpoint serviceEndpoint, System.ServiceModel.Description.ClientCredentials clientCredentials);
        
        public HeliocentricTrajectoriesInterfaceClient() : 
                base(HeliocentricTrajectoriesInterfaceClient.GetDefaultBinding(), HeliocentricTrajectoriesInterfaceClient.GetDefaultEndpointAddress())
        {
            this.Endpoint.Name = EndpointConfiguration.HeliocentricTrajectoriesPort.ToString();
            ConfigureEndpoint(this.Endpoint, this.ClientCredentials);
        }
        
        public HeliocentricTrajectoriesInterfaceClient(EndpointConfiguration endpointConfiguration) : 
                base(HeliocentricTrajectoriesInterfaceClient.GetBindingForEndpoint(endpointConfiguration), HeliocentricTrajectoriesInterfaceClient.GetEndpointAddress(endpointConfiguration))
        {
            this.Endpoint.Name = endpointConfiguration.ToString();
            ConfigureEndpoint(this.Endpoint, this.ClientCredentials);
        }
        
        public HeliocentricTrajectoriesInterfaceClient(EndpointConfiguration endpointConfiguration, string remoteAddress) : 
                base(HeliocentricTrajectoriesInterfaceClient.GetBindingForEndpoint(endpointConfiguration), new System.ServiceModel.EndpointAddress(remoteAddress))
        {
            this.Endpoint.Name = endpointConfiguration.ToString();
            ConfigureEndpoint(this.Endpoint, this.ClientCredentials);
        }
        
        public HeliocentricTrajectoriesInterfaceClient(EndpointConfiguration endpointConfiguration, System.ServiceModel.EndpointAddress remoteAddress) : 
                base(HeliocentricTrajectoriesInterfaceClient.GetBindingForEndpoint(endpointConfiguration), remoteAddress)
        {
            this.Endpoint.Name = endpointConfiguration.ToString();
            ConfigureEndpoint(this.Endpoint, this.ClientCredentials);
        }
        
        public HeliocentricTrajectoriesInterfaceClient(System.ServiceModel.Channels.Binding binding, System.ServiceModel.EndpointAddress remoteAddress) : 
                base(binding, remoteAddress)
        {
        }
        
        [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
        System.Threading.Tasks.Task<NASA.getAllObjectsResponse> NASA.HeliocentricTrajectoriesInterface.getAllObjectsAsync(NASA.getAllObjectsRequest request)
        {
            return base.Channel.getAllObjectsAsync(request);
        }
        
        public System.Threading.Tasks.Task<NASA.getAllObjectsResponse> getAllObjectsAsync()
        {
            NASA.getAllObjectsRequest inValue = new NASA.getAllObjectsRequest();
            return ((NASA.HeliocentricTrajectoriesInterface)(this)).getAllObjectsAsync(inValue);
        }
        
        [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
        System.Threading.Tasks.Task<NASA.getTrajectoriesResponse> NASA.HeliocentricTrajectoriesInterface.getTrajectoriesAsync(NASA.getTrajectoriesRequest request)
        {
            return base.Channel.getTrajectoriesAsync(request);
        }
        
        public System.Threading.Tasks.Task<NASA.getTrajectoriesResponse> getTrajectoriesAsync(NASA.request arg0)
        {
            NASA.getTrajectoriesRequest inValue = new NASA.getTrajectoriesRequest();
            inValue.arg0 = arg0;
            return ((NASA.HeliocentricTrajectoriesInterface)(this)).getTrajectoriesAsync(inValue);
        }
        
        [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
        System.Threading.Tasks.Task<NASA.getPrivacyAndImportantNoticesResponse> NASA.HeliocentricTrajectoriesInterface.getPrivacyAndImportantNoticesAsync(NASA.getPrivacyAndImportantNoticesRequest request)
        {
            return base.Channel.getPrivacyAndImportantNoticesAsync(request);
        }
        
        public System.Threading.Tasks.Task<NASA.getPrivacyAndImportantNoticesResponse> getPrivacyAndImportantNoticesAsync()
        {
            NASA.getPrivacyAndImportantNoticesRequest inValue = new NASA.getPrivacyAndImportantNoticesRequest();
            return ((NASA.HeliocentricTrajectoriesInterface)(this)).getPrivacyAndImportantNoticesAsync(inValue);
        }
        
        [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Advanced)]
        System.Threading.Tasks.Task<NASA.getAcknowledgementsResponse> NASA.HeliocentricTrajectoriesInterface.getAcknowledgementsAsync(NASA.getAcknowledgementsRequest request)
        {
            return base.Channel.getAcknowledgementsAsync(request);
        }
        
        public System.Threading.Tasks.Task<NASA.getAcknowledgementsResponse> getAcknowledgementsAsync()
        {
            NASA.getAcknowledgementsRequest inValue = new NASA.getAcknowledgementsRequest();
            return ((NASA.HeliocentricTrajectoriesInterface)(this)).getAcknowledgementsAsync(inValue);
        }
        
        public virtual System.Threading.Tasks.Task OpenAsync()
        {
            return System.Threading.Tasks.Task.Factory.FromAsync(((System.ServiceModel.ICommunicationObject)(this)).BeginOpen(null, null), new System.Action<System.IAsyncResult>(((System.ServiceModel.ICommunicationObject)(this)).EndOpen));
        }
        
        public virtual System.Threading.Tasks.Task CloseAsync()
        {
            return System.Threading.Tasks.Task.Factory.FromAsync(((System.ServiceModel.ICommunicationObject)(this)).BeginClose(null, null), new System.Action<System.IAsyncResult>(((System.ServiceModel.ICommunicationObject)(this)).EndClose));
        }
        
        private static System.ServiceModel.Channels.Binding GetBindingForEndpoint(EndpointConfiguration endpointConfiguration)
        {
            if ((endpointConfiguration == EndpointConfiguration.HeliocentricTrajectoriesPort))
            {
                System.ServiceModel.BasicHttpBinding result = new System.ServiceModel.BasicHttpBinding();
                result.MaxBufferSize = int.MaxValue;
                result.ReaderQuotas = System.Xml.XmlDictionaryReaderQuotas.Max;
                result.MaxReceivedMessageSize = int.MaxValue;
                result.AllowCookies = true;
                result.Security.Mode = System.ServiceModel.BasicHttpSecurityMode.Transport;
                return result;
            }
            throw new System.InvalidOperationException(string.Format("Не удалось найти конечную точку с именем \"{0}\".", endpointConfiguration));
        }
        
        private static System.ServiceModel.EndpointAddress GetEndpointAddress(EndpointConfiguration endpointConfiguration)
        {
            if ((endpointConfiguration == EndpointConfiguration.HeliocentricTrajectoriesPort))
            {
                return new System.ServiceModel.EndpointAddress("https://sscweb.gsfc.nasa.gov/WS/helio/1/HeliocentricTrajectoriesService");
            }
            throw new System.InvalidOperationException(string.Format("Не удалось найти конечную точку с именем \"{0}\".", endpointConfiguration));
        }
        
        private static System.ServiceModel.Channels.Binding GetDefaultBinding()
        {
            return HeliocentricTrajectoriesInterfaceClient.GetBindingForEndpoint(EndpointConfiguration.HeliocentricTrajectoriesPort);
        }
        
        private static System.ServiceModel.EndpointAddress GetDefaultEndpointAddress()
        {
            return HeliocentricTrajectoriesInterfaceClient.GetEndpointAddress(EndpointConfiguration.HeliocentricTrajectoriesPort);
        }
        
        public enum EndpointConfiguration
        {
            
            HeliocentricTrajectoriesPort,
        }
    }
}