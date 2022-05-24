﻿//------------------------------------------------------------------------------
// <auto-generated>
//     Этот код создан программой.
//
//     Изменения в этом файле могут привести к неправильной работе и будут потеряны в случае
//     повторного создания кода.
// </auto-generated>
//------------------------------------------------------------------------------

namespace CBR.Coins
{
    
    
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    [System.ServiceModel.ServiceContractAttribute(Namespace="http://web.cbr.ru/", ConfigurationName="CBR.Coins.CoinsBaseWSSoap")]
    public interface CoinsBaseWSSoap
    {
        
        [System.ServiceModel.OperationContractAttribute(Action="http://web.cbr.ru/EnumSeries", ReplyAction="*")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<CBR.Coins.ArrayOfXElement> EnumSeriesAsync();
        
        [System.ServiceModel.OperationContractAttribute(Action="http://web.cbr.ru/EnumSeriesXML", ReplyAction="*")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<System.Xml.XmlNode> EnumSeriesXMLAsync();
        
        [System.ServiceModel.OperationContractAttribute(Action="http://web.cbr.ru/EnumNominals", ReplyAction="*")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<CBR.Coins.ArrayOfXElement> EnumNominalsAsync();
        
        [System.ServiceModel.OperationContractAttribute(Action="http://web.cbr.ru/EnumNominalsXML", ReplyAction="*")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<System.Xml.XmlNode> EnumNominalsXMLAsync();
        
        [System.ServiceModel.OperationContractAttribute(Action="http://web.cbr.ru/EnumMetals", ReplyAction="*")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<CBR.Coins.ArrayOfXElement> EnumMetalsAsync();
        
        [System.ServiceModel.OperationContractAttribute(Action="http://web.cbr.ru/EnumMetalsXML", ReplyAction="*")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<System.Xml.XmlNode> EnumMetalsXMLAsync();
        
        [System.ServiceModel.OperationContractAttribute(Action="http://web.cbr.ru/SearchMonetXML", ReplyAction="*")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<System.Xml.XmlNode> SearchMonetXMLAsync(string SearchPhrase, int year, int nominal, int metal_id, int serie_id, int is_investment);
        
        [System.ServiceModel.OperationContractAttribute(Action="http://web.cbr.ru/SearchMonet", ReplyAction="*")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<CBR.Coins.ArrayOfXElement> SearchMonetAsync(string SearchPhrase, int year, int nominal, int metal_id, int serie_id, int is_investment);
        
        [System.ServiceModel.OperationContractAttribute(Action="http://web.cbr.ru/GetMonetDetailInfoXML", ReplyAction="*")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<System.Xml.XmlNode> GetMonetDetailInfoXMLAsync(string CatNumber, bool Eng);
        
        [System.ServiceModel.OperationContractAttribute(Action="http://web.cbr.ru/GetMonetDetailInfo", ReplyAction="*")]
        [System.ServiceModel.XmlSerializerFormatAttribute(SupportFaults=true)]
        System.Threading.Tasks.Task<CBR.Coins.ArrayOfXElement> GetMonetDetailInfoAsync(string CatNumber, bool Eng);
    }
    
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    public interface CoinsBaseWSSoapChannel : CBR.Coins.CoinsBaseWSSoap, System.ServiceModel.IClientChannel
    {
    }
    
    [System.Diagnostics.DebuggerStepThroughAttribute()]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("Microsoft.Tools.ServiceModel.Svcutil", "2.0.3")]
    public partial class CoinsBaseWSSoapClient : System.ServiceModel.ClientBase<CBR.Coins.CoinsBaseWSSoap>, CBR.Coins.CoinsBaseWSSoap
    {
        
        /// <summary>
        /// Реализуйте этот разделяемый метод для настройки конечной точки службы.
        /// </summary>
        /// <param name="serviceEndpoint">Настраиваемая конечная точка</param>
        /// <param name="clientCredentials">Учетные данные клиента.</param>
        static partial void ConfigureEndpoint(System.ServiceModel.Description.ServiceEndpoint serviceEndpoint, System.ServiceModel.Description.ClientCredentials clientCredentials);
        
        public CoinsBaseWSSoapClient(EndpointConfiguration endpointConfiguration) : 
                base(CoinsBaseWSSoapClient.GetBindingForEndpoint(endpointConfiguration), CoinsBaseWSSoapClient.GetEndpointAddress(endpointConfiguration))
        {
            this.Endpoint.Name = endpointConfiguration.ToString();
            ConfigureEndpoint(this.Endpoint, this.ClientCredentials);
        }
        
        public CoinsBaseWSSoapClient(EndpointConfiguration endpointConfiguration, string remoteAddress) : 
                base(CoinsBaseWSSoapClient.GetBindingForEndpoint(endpointConfiguration), new System.ServiceModel.EndpointAddress(remoteAddress))
        {
            this.Endpoint.Name = endpointConfiguration.ToString();
            ConfigureEndpoint(this.Endpoint, this.ClientCredentials);
        }
        
        public CoinsBaseWSSoapClient(EndpointConfiguration endpointConfiguration, System.ServiceModel.EndpointAddress remoteAddress) : 
                base(CoinsBaseWSSoapClient.GetBindingForEndpoint(endpointConfiguration), remoteAddress)
        {
            this.Endpoint.Name = endpointConfiguration.ToString();
            ConfigureEndpoint(this.Endpoint, this.ClientCredentials);
        }
        
        public CoinsBaseWSSoapClient(System.ServiceModel.Channels.Binding binding, System.ServiceModel.EndpointAddress remoteAddress) : 
                base(binding, remoteAddress)
        {
        }
        
        public System.Threading.Tasks.Task<CBR.Coins.ArrayOfXElement> EnumSeriesAsync()
        {
            return base.Channel.EnumSeriesAsync();
        }
        
        public System.Threading.Tasks.Task<System.Xml.XmlNode> EnumSeriesXMLAsync()
        {
            return base.Channel.EnumSeriesXMLAsync();
        }
        
        public System.Threading.Tasks.Task<CBR.Coins.ArrayOfXElement> EnumNominalsAsync()
        {
            return base.Channel.EnumNominalsAsync();
        }
        
        public System.Threading.Tasks.Task<System.Xml.XmlNode> EnumNominalsXMLAsync()
        {
            return base.Channel.EnumNominalsXMLAsync();
        }
        
        public System.Threading.Tasks.Task<CBR.Coins.ArrayOfXElement> EnumMetalsAsync()
        {
            return base.Channel.EnumMetalsAsync();
        }
        
        public System.Threading.Tasks.Task<System.Xml.XmlNode> EnumMetalsXMLAsync()
        {
            return base.Channel.EnumMetalsXMLAsync();
        }
        
        public System.Threading.Tasks.Task<System.Xml.XmlNode> SearchMonetXMLAsync(string SearchPhrase, int year, int nominal, int metal_id, int serie_id, int is_investment)
        {
            return base.Channel.SearchMonetXMLAsync(SearchPhrase, year, nominal, metal_id, serie_id, is_investment);
        }
        
        public System.Threading.Tasks.Task<CBR.Coins.ArrayOfXElement> SearchMonetAsync(string SearchPhrase, int year, int nominal, int metal_id, int serie_id, int is_investment)
        {
            return base.Channel.SearchMonetAsync(SearchPhrase, year, nominal, metal_id, serie_id, is_investment);
        }
        
        public System.Threading.Tasks.Task<System.Xml.XmlNode> GetMonetDetailInfoXMLAsync(string CatNumber, bool Eng)
        {
            return base.Channel.GetMonetDetailInfoXMLAsync(CatNumber, Eng);
        }
        
        public System.Threading.Tasks.Task<CBR.Coins.ArrayOfXElement> GetMonetDetailInfoAsync(string CatNumber, bool Eng)
        {
            return base.Channel.GetMonetDetailInfoAsync(CatNumber, Eng);
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
            if ((endpointConfiguration == EndpointConfiguration.CoinsBaseWSSoap))
            {
                System.ServiceModel.BasicHttpBinding result = new System.ServiceModel.BasicHttpBinding();
                result.MaxBufferSize = int.MaxValue;
                result.ReaderQuotas = System.Xml.XmlDictionaryReaderQuotas.Max;
                result.MaxReceivedMessageSize = int.MaxValue;
                result.AllowCookies = true;
                return result;
            }
            if ((endpointConfiguration == EndpointConfiguration.CoinsBaseWSSoap12))
            {
                System.ServiceModel.Channels.CustomBinding result = new System.ServiceModel.Channels.CustomBinding();
                System.ServiceModel.Channels.TextMessageEncodingBindingElement textBindingElement = new System.ServiceModel.Channels.TextMessageEncodingBindingElement();
                textBindingElement.MessageVersion = System.ServiceModel.Channels.MessageVersion.CreateVersion(System.ServiceModel.EnvelopeVersion.Soap12, System.ServiceModel.Channels.AddressingVersion.None);
                result.Elements.Add(textBindingElement);
                System.ServiceModel.Channels.HttpTransportBindingElement httpBindingElement = new System.ServiceModel.Channels.HttpTransportBindingElement();
                httpBindingElement.AllowCookies = true;
                httpBindingElement.MaxBufferSize = int.MaxValue;
                httpBindingElement.MaxReceivedMessageSize = int.MaxValue;
                result.Elements.Add(httpBindingElement);
                return result;
            }
            throw new System.InvalidOperationException(string.Format("Не удалось найти конечную точку с именем \"{0}\".", endpointConfiguration));
        }
        
        private static System.ServiceModel.EndpointAddress GetEndpointAddress(EndpointConfiguration endpointConfiguration)
        {
            if ((endpointConfiguration == EndpointConfiguration.CoinsBaseWSSoap))
            {
                return new System.ServiceModel.EndpointAddress("http://www.cbr.ru/CoinsBaseWS/CoinsBaseWS.asmx");
            }
            if ((endpointConfiguration == EndpointConfiguration.CoinsBaseWSSoap12))
            {
                return new System.ServiceModel.EndpointAddress("http://www.cbr.ru/CoinsBaseWS/CoinsBaseWS.asmx");
            }
            throw new System.InvalidOperationException(string.Format("Не удалось найти конечную точку с именем \"{0}\".", endpointConfiguration));
        }
        
        public enum EndpointConfiguration
        {
            
            CoinsBaseWSSoap,
            
            CoinsBaseWSSoap12,
        }
    }
    
    [System.Xml.Serialization.XmlSchemaProviderAttribute(null, IsAny=true)]
    [System.CodeDom.Compiler.GeneratedCodeAttribute("dotnet-svcutil-lib", "2.0.3.0")]
    public partial class ArrayOfXElement : object, System.Xml.Serialization.IXmlSerializable
    {
        
        private System.Collections.Generic.List<System.Xml.Linq.XElement> nodesList = new System.Collections.Generic.List<System.Xml.Linq.XElement>();
        
        public ArrayOfXElement()
        {
        }
        
        public virtual System.Collections.Generic.List<System.Xml.Linq.XElement> Nodes
        {
            get
            {
                return this.nodesList;
            }
        }
        
        public virtual System.Xml.Schema.XmlSchema GetSchema()
        {
            throw new System.NotImplementedException();
        }
        
        public virtual void WriteXml(System.Xml.XmlWriter writer)
        {
            System.Collections.Generic.IEnumerator<System.Xml.Linq.XElement> e = nodesList.GetEnumerator();
            for (
            ; e.MoveNext(); 
            )
            {
                ((System.Xml.Serialization.IXmlSerializable)(e.Current)).WriteXml(writer);
            }
        }
        
        public virtual void ReadXml(System.Xml.XmlReader reader)
        {
            for (
            ; (reader.NodeType != System.Xml.XmlNodeType.EndElement); 
            )
            {
                if ((reader.NodeType == System.Xml.XmlNodeType.Element))
                {
                    System.Xml.Linq.XElement elem = new System.Xml.Linq.XElement("default");
                    ((System.Xml.Serialization.IXmlSerializable)(elem)).ReadXml(reader);
                    Nodes.Add(elem);
                }
                else
                {
                    reader.Skip();
                }
            }
        }
    }
}